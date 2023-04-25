/*
  Input:
  Global: radious, similarity, sortedStrings.
  Output:
  Global: searchResults.
*/
namespace layered 
{
    int layeredMaxMoreThresholds[MAXIMUM_K_MER_LENGTH];
    int layeredMinMoreThresholds[MAXIMUM_K_MER_LENGTH];

    int **allQueryCounts = 0;

    int layeredClusterSearch(Index index, int k, IndexVector &clusterSearchResults)
    {
    
        int maxMore = 0;
        int minMore = 0;
        int numberOfKmers = cluster_information::numberOfKmers[k];
        KmerCount *minCounts = &cluster_information::allMinCounts[k][index * numberOfKmers];
        KmerCount *maxCounts = &cluster_information::allMaxCounts[k][index * numberOfKmers];
        int minMoreThreshold = layeredMinMoreThresholds[k];
        int maxMoreThreshold = layeredMaxMoreThresholds[k];
        int *queryCounts = allQueryCounts[k];
        {
            int minCount;
            int maxCount;
            int queryCount;
            for (Index i = 0; i < numberOfKmers; ++i) 
            {
                maxCount = maxCounts[i];
                queryCount = queryCounts[i];
                if (maxCount > queryCount) 
                {
                    maxMore += maxCount - queryCount;
                    minCount = minCounts[i];
                    if (minCount > queryCount) 
                    { 
                        minMore += minCount - queryCount;
                        if (minMore > minMoreThreshold) {
                            return 0;
                        } 
                        if (approximateFilter2 && minMore * (cluster_information::r[index] - cluster_information::l[index]) > minMoreThreshold)
                            return 0;
                    } 
                } 
            } // for (Index i = 0; i < numberOfKmers; ++i)
        }
        if (maxMore <= maxMoreThreshold) 
        {
            if (k + 1 < k_mer_length) 
            {
                return layeredClusterSearch(index, k + 1, clusterSearchResults);
            }
            for (Index i = cluster_information::l[index]; i < cluster_information::r[index]; ++i)
                clusterSearchResults.push_back(i);
            return cluster_information::r[index] - cluster_information::l[index];
        }
        return (layeredClusterSearch(index + 1, k, clusterSearchResults) + layeredClusterSearch(index + 2 * cluster_information::sizeOfLeftCluster[index], k, clusterSearchResults));
    }
} 

void search_without_backpointer(const PositionType l, const PositionType r, PositionType pos, const Index intervalIndex, const uint16_t tableIndex, const int globalOffset)
{
    DpTable *table = tables[tableIndex];
    std::vector<CharType> *localWord = &word[tableIndex];
    std::vector<PositionType> *localWordMinLength = &wordMinLength[tableIndex];
    std::vector<PositionType> *localLeft = &left[tableIndex];
    std::vector<PositionType> *localRight = &right[tableIndex];
    std::vector<PositionType> *localMinCost = &minCost[tableIndex];

    CostType maxLengthRadious = std::min(static_cast<CostType>(floor((1 - similarity) * sortedStringsLengthsMaxTree[tableIndex][intervalIndex])), radious);
    if (mismatches >= 0)
      maxLengthRadious = static_cast<CostType>(radious);

    PositionType minSortedStringsLength = sortedStringsLengthsMinTree[tableIndex][intervalIndex];
    PositionType maxRight = querySequenceLength - minSortedStringsLength + radious;

    do 
    {
        if (sortedStrings[globalOffset + l][pos] != sortedStrings[globalOffset + r][pos]) 
        {
            PositionType mid = (l + r) / 2;
            Index leftIntervalIndex = intervalIndex + 1;
            Index rightIntervalIndex = 2 * (mid - l + 1) + intervalIndex;
            if (sortedStringsLengthsMinTree[tableIndex][leftIntervalIndex] < sortedStringsLengthsMinTree[tableIndex][rightIntervalIndex]) 
            {
              search_without_backpointer(l,       mid, pos, leftIntervalIndex, tableIndex, globalOffset);
              search_without_backpointer(mid + 1, r,   pos, rightIntervalIndex, tableIndex, globalOffset);
            } 
            else 
            {
              search_without_backpointer(mid + 1, r,   pos, rightIntervalIndex, tableIndex, globalOffset);
              search_without_backpointer(l,       mid, pos, leftIntervalIndex, tableIndex, globalOffset);
            }
            return;
        } 
        if (endOfStringChar == sortedStrings[globalOffset + l][pos]) 
        {
            SearchResult result;
            PositionType minPos = std::max(0, (*localLeft)[pos] + pos);
            for(PositionType i = minPos + 1; i <= (*localRight)[pos] + pos; ++i)
                if ((*table)(pos, i) < (*table)(pos, minPos))
                    minPos = i;
          
            // Check the number of differences is less than or equal to similarity times the length of the shorter sequence.
            float acceptable_differences = pos * (1 - similarity);
            if (mismatches >= 0)
                acceptable_differences = mismatches;
            if ((*table)(pos, minPos) <= acceptable_differences) 
            {
                for (PositionType i = l; i <= r; ++i) 
                {
                    result.number = i + globalOffset;
                    threadSearchResults[tableIndex].push_back(result);
                } 
            } 
            return;
        }
    
        if ((*localWord)[pos + 1] != sortedStrings[globalOffset + l][pos] || (*localWordMinLength)[pos + 1] > minSortedStringsLength) 
        {
            (*localWord)[pos + 1] = sortedStrings[globalOffset + l][pos];
            (*localWord)[pos + 2] = invalidChar;
            (*localWordMinLength)[pos + 1] = minSortedStringsLength;
            {  // Update the row table(pos + 1, *)
                const PositionType minColumn = std::max(0, (*localLeft)[pos] + (pos + 1));
                const PositionType maxColumn = std::min(querySequenceLength, (*localRight)[pos] + (pos + 1));
                CostType currentMinCost = MAX_COST;
                CostType currentCost;
                const PositionType i = pos + 1;
                const CharType c1 = (*localWord)[i];
                CharType c2; 
                CostType matchCost, gapCost2, gapCost1;
                CostType leftCost, downLeftCost, downCost;
                leftCost = MAX_COST;
                if (minColumn > 0)
                    downLeftCost = (*table)(i - 1, minColumn - 1);
                else 
                    downLeftCost = MAX_COST;
                for (PositionType j = minColumn; j <= maxColumn; ++j) 
                {
                    { // Update table(i, j) = table(pos + 1, j)
                        currentCost = MAX_COST;
                        c2 = querySequencePointer[j - 1];
                        matchCost = downLeftCost + cost(c1, c2);
                        if (matchCost < currentCost) 
                            currentCost = matchCost;
                        gapCost2 = leftCost + cost(gapChar, c2);
                        if (gapCost2 < currentCost) 
                            currentCost = gapCost2;
                        downCost = (*table)(i - 1, j);
                        gapCost1 = downCost + cost(c1, gapChar);
                        if (gapCost1 < currentCost) 
                            currentCost = gapCost1;
                        (*table)(i, j) = currentCost;
                        leftCost = currentCost;
                        downLeftCost = downCost;
                    } // Update table(i, j) = table(pos + 1, j)
                    if (currentCost < currentMinCost)
                        currentMinCost = currentCost;
                } // for (PositionType j...
                (*localMinCost)[pos + 1] = currentMinCost;
            }
            if ((*localRight)[pos] + (pos + 1) < querySequenceLength)
                (*table)(pos + 1, (*localRight)[pos] + (pos + 1) + 1) = MAX_COST;
            { // Find left[pos + 1].
                PositionType j = std::max(-(pos + 1), (*localLeft)[pos]);
                PositionType maxIndex = std::min(querySequenceLength - (pos + 1), (*localRight)[pos]);
                while ((j <= maxIndex) && ((*table)(pos + 1, j + (pos + 1)) > radious))
                    ++j;
                if (j == -(pos + 1))
                    (*localLeft)[pos + 1] = j - 1;
                else
                    (*localLeft)[pos + 1] = j;
            }
            { // Find right[pos + 1].
                PositionType j = std::min(querySequenceLength - (pos + 1), (*localRight)[pos]);
                PositionType minIndex = std::max(-(pos + 1), (*localLeft)[pos]);
                while ((j >= minIndex) && ((*table)(pos + 1, j + (pos + 1)) > radious))
                    --j;
                (*localRight)[pos + 1] = std::min(j, maxRight);
            } // Find right[pos + 1]
        } // if (word[pos + 1] != sortedStrings[l][pos])
    
    }while ((*localMinCost)[++pos] <= maxLengthRadious);
}

void *threaded_search_without_backpointer(void *read_args) 
{
    struct thread_args *args = static_cast<struct thread_args *>(read_args);
    const PositionType l = args->l;
    const PositionType r = args->r;
    PositionType pos = args->pos;
    const Index intervalIndex = args->intervalIndex;
    const uint16_t tableIndex = args->tableIndex;
    const int globalOffset = args->globalOffset;
    search_without_backpointer(l, r, pos, intervalIndex, tableIndex, globalOffset);
    return NULL;
}