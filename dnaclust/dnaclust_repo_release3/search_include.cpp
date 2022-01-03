
#undef SEARCH_FUNCTION_NAME

#if ALIGN_WITH_BACKPOINTER
#define SEARCH_FUNCTION_NAME search
#else // ALIGN_WITH_BACKPOINTER is false
#define SEARCH_FUNCTION_NAME search_without_backpointer
#endif // ALIGN_WITH_BACKPOINTER

/*
  Input:
  Global: radious, similarity, sortedStrings.

  Output:
  Global: searchResults.
*/


void SEARCH_FUNCTION_NAME(const PositionType l, const PositionType r, PositionType pos, const Index intervalIndex, const uint16_t tableIndex, const int globalOffset)
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

  do {


    if (sortedStrings[globalOffset + l][pos] != sortedStrings[globalOffset + r][pos]) {
      PositionType mid = (l + r) / 2;
      Index leftIntervalIndex = intervalIndex + 1;
      Index rightIntervalIndex = 2 * (mid - l + 1) + intervalIndex;
      if (sortedStringsLengthsMinTree[tableIndex][leftIntervalIndex] < sortedStringsLengthsMinTree[tableIndex][rightIntervalIndex]) {
	SEARCH_FUNCTION_NAME(l,       mid, pos, leftIntervalIndex, tableIndex, globalOffset);
	SEARCH_FUNCTION_NAME(mid + 1, r,   pos, rightIntervalIndex, tableIndex, globalOffset);
      } else {
	SEARCH_FUNCTION_NAME(mid + 1, r,   pos, rightIntervalIndex, tableIndex, globalOffset);
	SEARCH_FUNCTION_NAME(l,       mid, pos, leftIntervalIndex, tableIndex, globalOffset);
      }
      return;
    } 

    if (endOfStringChar == sortedStrings[globalOffset + l][pos]) {
      SearchResult result;

      PositionType minPos = std::max(0, (*localLeft)[pos] + pos);
      for(PositionType i = minPos + 1; i <= (*localRight)[pos] + pos; ++i)
	if ((*table)(pos, i) < (*table)(pos, minPos))
	  minPos = i;
      
      // Check the number of differences is less than or equal to similarity times the length of the shorter sequence.
      float acceptable_differences = pos * (1 - similarity);
      if (mismatches >= 0)
        acceptable_differences = mismatches;
      if ((*table)(pos, minPos) <= acceptable_differences) {

	if (outputMultipleAlignment) {
	  boost::tie(result.cost, result.func) = findAlignmentFunction(pos, querySequenceLength, minPos);
	}

	// Try to lower the amount of time we try to claim the lock.
	//if (l <= r){
		//pthread_mutex_lock(&running_mutex);
		for (PositionType i = l; i <= r; ++i) {
		  result.number = i + globalOffset;
		  threadSearchResults[tableIndex].push_back(result);
		  //++resultsCounter[tableIndex];
		} // for (PositionType i = l;; i <= r; ++i)

      } // if (table(pos, minPos) <= pos * (1 - similarity))

      return;
    }
    

    if ((*localWord)[pos + 1] != sortedStrings[globalOffset + l][pos]
	|| (*localWordMinLength)[pos + 1] > minSortedStringsLength
	) {
      // Row table(pos + 1, *) must be updated.
      // Thus minCost[pos + 1] must be updated.

      (*localWord)[pos + 1] = sortedStrings[globalOffset + l][pos];
      (*localWord)[pos + 2] = invalidChar;

      (*localWordMinLength)[pos + 1] = minSortedStringsLength;

      {  // Update the row table(pos + 1, *)
	const PositionType minColumn = std::max(0, (*localLeft)[pos] + (pos + 1));
	const PositionType maxColumn = std::min(querySequenceLength, (*localRight)[pos] + (pos + 1));
	CostType currentMinCost = MAX_COST;


	CostType currentCost;
#if ALIGN_WITH_BACKPOINTER
	BackPointer currentPointer;
#endif // ALIGN_WITH_BACKPOINTER
	const PositionType i = pos + 1;
	const CharType c1 = (*localWord)[i];
	CharType c2; // = querySequencePointer[j - 1]
	CostType matchCost, gapCost2, gapCost1;
	
	CostType leftCost, downLeftCost, downCost;


	leftCost = MAX_COST;
	if (minColumn > 0) {
	  downLeftCost = (*table)(i - 1, minColumn - 1);
	} else {
	  downLeftCost = MAX_COST;
	}
	

	for (PositionType j = minColumn; j <= maxColumn; ++j) {

	  { // Update table(i, j) = table(pos + 1, j)
	    
#if ALIGN_WITH_BACKPOINTER
	    currentPointer = ROOT;
#endif // ALIGN_WITH_BACKPOINTER
	    currentCost = MAX_COST;
	    
	    // c1 == word[i]
	    c2 = querySequencePointer[j - 1];
	    
	    matchCost = downLeftCost + cost(c1, c2);
	    if (matchCost < currentCost) {
	      currentCost = matchCost;
#if ALIGN_WITH_BACKPOINTER
	      currentPointer = DOWN_LEFT;
#endif // ALIGN_WITH_BACKPOINTER
	    }
	      
	    gapCost2 = leftCost + cost(gapChar, c2);
	    if (gapCost2 < currentCost) {
	      currentCost = gapCost2;
#if ALIGN_WITH_BACKPOINTER
	      currentPointer = LEFT;
#endif // ALIGN_WITH_BACKPOINTER
	    }
	    
	    downCost = (*table)(i - 1, j);

	    gapCost1 = downCost + cost(c1, gapChar);
	    if (gapCost1 < currentCost) {
	      currentCost = gapCost1;
#if ALIGN_WITH_BACKPOINTER
	      currentPointer = DOWN;
#endif // ALIGN_WITH_BACKPOINTER
	    }
	    
	    (*table)(i, j) = currentCost;
#if ALIGN_WITH_BACKPOINTER
	    backTable(i, j) = currentPointer;
#endif // ALIGN_WITH_BACKPOINTER

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
    
  } while ((*localMinCost)[++pos] <= maxLengthRadious);

}
