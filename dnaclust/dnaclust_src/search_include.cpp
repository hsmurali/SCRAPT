
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


void SEARCH_FUNCTION_NAME(const PositionType l, const PositionType r, PositionType pos, const Index intervalIndex)
{

  CostType maxLengthRadious = std::min(static_cast<CostType>(floor((1 - similarity) * sortedStringsLengthsMaxTree[intervalIndex])), radious);
  PositionType minSortedStringsLength = sortedStringsLengthsMinTree[intervalIndex];
  PositionType maxRight = querySequenceLength - minSortedStringsLength + radious;

  do {


    if (sortedStrings[l][pos] != sortedStrings[r][pos]) {
      PositionType mid = (l + r) / 2;
      Index leftIntervalIndex = intervalIndex + 1;
      Index rightIntervalIndex = 2 * (mid - l + 1) + intervalIndex;
      if (sortedStringsLengthsMinTree[leftIntervalIndex] < sortedStringsLengthsMinTree[rightIntervalIndex]) {
	SEARCH_FUNCTION_NAME(l,       mid, pos, leftIntervalIndex);
	SEARCH_FUNCTION_NAME(mid + 1, r,   pos, rightIntervalIndex);
      } else {
	SEARCH_FUNCTION_NAME(mid + 1, r,   pos, rightIntervalIndex);
	SEARCH_FUNCTION_NAME(l,       mid, pos, leftIntervalIndex);
      }
      return;
    } 

    if (endOfStringChar == sortedStrings[l][pos]) {
      SearchResult result;

      PositionType minPos = std::max(0, left[pos] + pos);
      for(PositionType i = minPos + 1; i <= right[pos] + pos; ++i)
	if (table(pos, i) < table(pos, minPos))
	  minPos = i;
      
      // Check the number of differences is less than or equal to similarity times the length of the shorter sequence.
      if (table(pos, minPos) <= pos * (1 - similarity)) {

	if (outputMultipleAlignment) {
	  boost::tie(result.cost, result.func) = findAlignmentFunction(pos, querySequenceLength, minPos);
	}

	for (PositionType i = l; i <= r; ++i) {
	  result.number = i;
	  searchResults.push_back(result);
	} // for (PositionType i = l;; i <= r; ++i)

      } // if (table(pos, minPos) <= pos * (1 - similarity))

      return;
    }
    

    if (word[pos + 1] != sortedStrings[l][pos]
	|| wordMinLength[pos + 1] > minSortedStringsLength
	) {
      // Row table(pos + 1, *) must be updated.
      // Thus minCost[pos + 1] must be updated.


      word[pos + 1] = sortedStrings[l][pos];
      word[pos + 2] = invalidChar;

      wordMinLength[pos + 1] = minSortedStringsLength;

      {  // Update the row table(pos + 1, *)
	const PositionType minColumn = std::max(0, left[pos] + (pos + 1));
	const PositionType maxColumn = std::min(querySequenceLength, right[pos] + (pos + 1));
	CostType currentMinCost = MAX_COST;


	CostType currentCost;
#if ALIGN_WITH_BACKPOINTER
	BackPointer currentPointer;
#endif // ALIGN_WITH_BACKPOINTER
	const PositionType i = pos + 1;
	const CharType c1 = word[i];
	CharType c2; // = querySequencePointer[j - 1]
	CostType matchCost, gapCost2, gapCost1;
	
	CostType leftCost, downLeftCost, downCost;


	leftCost = MAX_COST;
	if (minColumn > 0) {
	  downLeftCost = table(i - 1, minColumn - 1);
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
	    
	    downCost = table(i - 1, j);

	    gapCost1 = downCost + cost(c1, gapChar);
	    if (gapCost1 < currentCost) {
	      currentCost = gapCost1;
#if ALIGN_WITH_BACKPOINTER
	      currentPointer = DOWN;
#endif // ALIGN_WITH_BACKPOINTER
	    }
	    
	    table(i, j) = currentCost;
#if ALIGN_WITH_BACKPOINTER
	    backTable(i, j) = currentPointer;
#endif // ALIGN_WITH_BACKPOINTER

	    leftCost = currentCost;
	    downLeftCost = downCost;
	    
	  } // Update table(i, j) = table(pos + 1, j)

	  
	  if (currentCost < currentMinCost)
	    currentMinCost = currentCost;
	} // for (PositionType j...


	minCost[pos + 1] = currentMinCost;
      }

      if (right[pos] + (pos + 1) < querySequenceLength)
	table(pos + 1, right[pos] + (pos + 1) + 1) = MAX_COST;
      
      
      { // Find left[pos + 1].
	PositionType j = std::max(-(pos + 1), left[pos]);
	PositionType maxIndex = std::min(querySequenceLength - (pos + 1), right[pos]);
	while ((j <= maxIndex) && (table(pos + 1, j + (pos + 1)) > radious))
	  ++j;
	if (j == -(pos + 1))
	  left[pos + 1] = j - 1;
	else
	  left[pos + 1] = j;
      }

      { // Find right[pos + 1].
	PositionType j = std::min(querySequenceLength - (pos + 1), right[pos]);
	PositionType minIndex = std::max(-(pos + 1), left[pos]);
	while ((j >= minIndex) && (table(pos + 1, j + (pos + 1)) > radious))
	  --j;
	
	right[pos + 1] = std::min(j, maxRight);

      } // Find right[pos + 1]
    
    } // if (word[pos + 1] != sortedStrings[l][pos])
    
  } while (minCost[++pos] <= maxLengthRadious);

}
