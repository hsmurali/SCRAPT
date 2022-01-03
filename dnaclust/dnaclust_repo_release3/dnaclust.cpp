#include <sstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <stack>
#include <fstream>
// #include <boost/bind.hpp>
// #include <boost/date_time/posix_time/posix_time.hpp>
// #include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/program_options.hpp>		
#include <stdint.h>
//#include <boost/thread.hpp>
//#include <boost/thread/condition_variable.hpp>
#include <pthread.h>
namespace po = boost::program_options;

#include "ternary_sort.hpp"
#include "fasta.hpp"

#define MATRIX_RESIZABLE 1
#include "multi_dim.hpp"
#include "utility.hpp"


typedef char CharType;
typedef const CharType * ConstCharPointer;
typedef std::vector<ConstCharPointer> ConstCharPointerVector;
typedef int PositionType;
typedef int SequenceNumber;
typedef intptr_t Index;

uint16_t running_threads = 0;
size_t num_threads = 1;
int threads_completed = 0;
pthread_mutex_t running_mutex;

struct thread_args
 {
    PositionType l;
    PositionType r;
    PositionType pos;
    Index intervalIndex;
    uint16_t tableIndex;
    int globalOffset;
};

SequenceNumber numberOfInputSequences;

bool outputMultipleAlignment = false;
bool leftGapsAllowed = false;
bool noKmerFilter = true;
bool kmerFilter = false;
#if OLD_CODE
bool approximateFilter = false;
double approximationFactor;
#endif // OLD_CODE
bool approximateFilter2 = false;
bool noOverlap = false;
bool recruit_only = false;
bool useFullQueryHeader = false;

ConstCharPointerVector sortedStrings;

ConstCharPointerVector lexicographicSequencePointers;
ConstCharPointerVector headerPointers;
typedef std::vector<PositionType> PositionVector;
PositionVector headerLengths;
PositionVector sequenceLengths;


ConstCharPointer querySequencePointer;
PositionType querySequenceLength;
ConstCharPointer queryHeaderPointer;
PositionType queryHeaderLength;


const CharType endOfStringChar = '\0';
const CharType invalidChar = '^';
const CharType gapChar = '-';

const char HEADER_INDICATOR = '%';
const char CLUSTER_INDICATOR = '#';

typedef int32_t CostType;
const CostType MAX_COST = std::numeric_limits<CostType>::max() / 2;

struct Context {
  int id;
};

const size_t MAX_THREADS = 60;
pthread_t searchWorkers[MAX_THREADS];
struct Context searchWorkersContext[MAX_THREADS];
struct thread_args searchWorkersArgs[MAX_THREADS];
pthread_mutex_t searchWorkersMutexes[MAX_THREADS];
pthread_cond_t searchWorkersConds[MAX_THREADS];
pthread_cond_t searchWorkersCompleted;
pthread_mutex_t allThreadsMutex;

const size_t MAX_LENGTH = 4500;
std::vector< std::vector<CharType> > word(MAX_THREADS, std::vector<CharType>(MAX_LENGTH + 1));
std::vector< std::vector<PositionType> > wordMinLength(MAX_THREADS, std::vector<PositionType>(MAX_LENGTH + 1));
std::vector< std::vector<CostType> > minCost(MAX_THREADS, std::vector<CostType>(MAX_LENGTH + 1));
std::vector< std::vector<PositionType> > left(MAX_THREADS, std::vector<PositionType>(MAX_LENGTH + 1));
std::vector< std::vector<PositionType> > right(MAX_THREADS, std::vector<PositionType>(MAX_LENGTH + 1));

CostType radious;

float DEFAULT_SIMILARITY = 0.99;
float similarity;

float DEFAULT_MISMATCHES = -1.0;
float mismatches;

bool printInvertedIndex = false;
bool assignAmbiguous = false;
std::map< std::string, std::vector<std::string> > final_clusters;
std::map< std::string, std::vector<std::string> > inverted_index;
std::map< std::string, int > center_nonambig_counts;

template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"{ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, " "));
    os<<"}";
    return os;
}

template <typename BaseType, typename ExponentType>
BaseType power(const BaseType base, const ExponentType exponent)
{
  BaseType result = 1;
  for (ExponentType i = 0; i < exponent; ++i)
    result *= base;
  return result;
}


const int MAXIMUM_K_MER_LENGTH = 6;
#if OLD_CODE
const int DEFAULT_K_MER_LENGTH = 4;
#endif // OLD_CODE
int k_mer_length;
int random_seed = 0;

enum Nucleotides{A, C, G, T, NUMBER_OF_NUCLEOTIDES};

typedef int16_t KmerNumber;

#if FIXED_K
const KmerNumber NUMBER_OF_K_MERS = power(static_cast<int>(NUMBER_OF_NUCLEOTIDES), K_MER_LENGTH);
#endif // FIXED_K
KmerNumber number_of_k_mers;

// Exception thrown when program gets an unkown nucleotide letter.
class UnknownNucleotide : public std::exception
{
public:
  UnknownNucleotide(CharType nucleotide) throw()
    : nucleotide_(nucleotide)
  { }

  virtual const char *what() const throw()
  {
    std::ostringstream messageStream;
    messageStream << "Sequence contains unknown nucleotide: " << nucleotide_ << ".";
    return messageStream.str().c_str();
  }

private:
  const CharType nucleotide_;
};

// Return an integer number in range 0..(NUMBER_OF_NUCLEOTIDES - 1) for a given nucleotide character.
int numberFromNucleotide(const CharType nucleotide)
{
  switch(tolower(nucleotide)) {
  case 'a':
    return A;
  case 'c':
    return C;
  case 'g':
    return G;
  case 't':
    return T;
  default:
    return A;
    throw UnknownNucleotide(nucleotide);
  }
}

typedef uint16_t KmerCount;
const KmerCount MAX_K_MER_COUNT = std::numeric_limits<KmerCount>::max();

// k-mer spectrum is a vector of integers counting how many of each k-mer there are in a sequence.
typedef std::vector<KmerCount> KmerSpectrum;

// Count the number of each k-mer in a sequence and return a vector of counts.
KmerSpectrum countKmers(ConstCharPointer sequence, PositionType length)
{
  KmerSpectrum spectrum(number_of_k_mers);

  if (length >= k_mer_length) {
    int i = 0;
    int kMerNumber = 0;
    for (; i < k_mer_length; ++i) {
      kMerNumber *= NUMBER_OF_NUCLEOTIDES;
      kMerNumber += numberFromNucleotide(sequence[i]);
    }
    
    ++spectrum[kMerNumber];
    
    for (; i < length; ++i) {
      kMerNumber *= NUMBER_OF_NUCLEOTIDES;
      kMerNumber %= number_of_k_mers;
      kMerNumber += numberFromNucleotide(sequence[i]);
      
      if (spectrum[kMerNumber] < MAX_K_MER_COUNT)
	++spectrum[kMerNumber];
    }
      
    
  }
  
  return spectrum;
}



typedef std::vector<Index> IndexVector;

IndexVector lexicographicIndexOfSortedStringsIndexes;

typedef multi_dim::Matrix<KmerCount> CountMatrix;
CountMatrix spectrumMatrix;

IndexVector lexicographicIndexOfSpectrumIndexes;


IndexVector spectrumSearchResults;

KmerSpectrum querySpectrum;
KmerSpectrum queryLowSpectrum;

IndexVector spectrumSplitGuide;
PositionVector spectrumRemainingGuide;

/*
  Search for potential matches based on their k-mer spectrum.
  The query's k-mer spectrum should in global variable 'querySpectrum'.
  The k-mer spectrums for the sequences to be clustered are in global variable 'specturms'.
  Results will be in global variable ''.
*/




typedef std::pair<KmerSpectrum, Index> SpectrumIndexPair;
typedef std::vector<SpectrumIndexPair> SpectrumIndexVector;
    
SpectrumIndexVector spectrumIndexes;


enum BackPointer {LEFT, DOWN, DOWN_LEFT, ROOT};

typedef multi_dim::Matrix<CostType> DpTable;

DpTable *tables[MAX_THREADS];
//DpTable table(MAX_LENGTH + 1, MAX_LENGTH + 1);

typedef multi_dim::Matrix<BackPointer> BackPointerTable;

BackPointerTable backTable(MAX_LENGTH + 1, MAX_LENGTH + 1);


inline const CostType cost(const CharType c1, const CharType c2)
{
  return (c1 != c2) ? 1 : 0;
}

#if COST_MATRIX
namespace cost_matrix {
  const size_t MAX_CHAR = 256;
  typedef multi_dim::Array2D<CostType, MAX_CHAR, MAX_CHAR> CostMatrix;

  CostMatrix cost;
  
}
#endif // COST_MATRIX

void initialize()
{
  // initialize() implementation depends on whether we allow gaps on both ends.
  
  for(uint16_t i = 0; i < num_threads; ++i){
    if (mismatches >= 0)
      left[i][0] = - static_cast<PositionType>(mismatches);
    else
      left[i][0] = - static_cast<PositionType>(floor(MAX_LENGTH * (1 - similarity)));
    right[i][0] = std::numeric_limits<PositionType>::max() / 2;
    std::fill(word[i].begin(), word[i].end(), invalidChar);

    tables[i] = new DpTable(MAX_LENGTH + 1, MAX_LENGTH + 1);
    (*tables[i])(0, 0) = 0;
    backTable(0, 0) = ROOT;

    for (size_t j = 1; j <  MAX_LENGTH + 1; ++j) {
      (*tables[i])(0, j) = leftGapsAllowed ? 0 : (*tables[i])(0, j - 1) + cost(gapChar, invalidChar);    
      backTable(0, j) = LEFT;
    }
  }
} // void initialize()


typedef std::vector<PositionType> AlignmentFunction;
typedef boost::tuple<CostType, AlignmentFunction> AlignmentCostAndFunction;


AlignmentCostAndFunction findAlignmentFunction(const PositionType len1, const PositionType len2, const PositionType pos2)
{

  AlignmentFunction alignmentFunc(len2 + 1);

  // TODO(cmhill): Remove
  CostType minCost = 0; // table(len1, pos2);


  PositionType i = len1;
  PositionType j = len2;

  assert(minCost <= radious);

  while (j > pos2){
    alignmentFunc[j] = i;
    --j;
  }

  while (backTable(i, j) != ROOT){

    switch(backTable(i, j)){
    case DOWN:
      alignmentFunc[j] = i;
      --i;
      break;
    case LEFT:
      alignmentFunc[j] = i;
      --j;
      break;
    case DOWN_LEFT:
      alignmentFunc[j] = i;
      --i;
      --j;
      break;
    case ROOT:
    default:
      std::stringstream errorStream;
      errorStream << "Invalid back pointer value; backTable(" << i << "," << j << ")=" << backTable(i, j) << ".";
      throw std::logic_error(errorStream.str());
    } // switch
  } // while (backTable(i, j) != ROOT)

  alignmentFunc[j] = i;
    
  return boost::make_tuple(minCost, alignmentFunc);
} // findAlignmentFunction2

struct SearchResult
{
  SequenceNumber number;
  CostType cost;
  AlignmentFunction func;
};

typedef std::vector<SearchResult> SearchResultVector;

SearchResultVector searchResults;

std::vector< SearchResultVector > threadSearchResults(MAX_THREADS);
std::vector< PositionVector > sortedStringsLengthsMinTree(MAX_THREADS);
std::vector< PositionVector > sortedStringsLengthsMaxTree(MAX_THREADS);
int resultsCounter[MAX_THREADS];

template <class Iterator1, class Iterator2>
Iterator2 buildMaxTree(Iterator1 first, Iterator1 last, Iterator2 intervalIterator)
{
  if (last > first) {
    if (last == first + 1) {
      *intervalIterator = *first;
    } else {
      Iterator1 mid = first + (last - 1 - first) / 2 + 1;
      *intervalIterator = std::max(*buildMaxTree(first,  mid, intervalIterator + 1),
				   *buildMaxTree(mid,   last, intervalIterator + 2 * (mid - first)));

    }
  }
  return intervalIterator;
}


template <class Iterator1, class Iterator2>
Iterator2 buildMinTree(Iterator1 first, Iterator1 last, Iterator2 intervalIterator)
{
  if (last > first) {
    if (last == first + 1) {
      *intervalIterator = *first;
    } else {
      Iterator1 mid = first + (last - 1 - first) / 2 + 1;
      *intervalIterator = std::min(*buildMinTree(first,  mid, intervalIterator + 1),
				   *buildMinTree(mid,   last, intervalIterator + 2 * (mid - first)));

    }
  }
  return intervalIterator;
}


/*
  The implementation of search function is seperated in another file.
  The reason for this is to reuse the code to make two search functions:
  The first search function does not keep track of back pointers and would not find the alignment but is supposedly faster.
  The second search function will find the alignment.
*/

#define ALIGN_WITH_BACKPOINTER 0
#include "search_include.cpp"

#undef ALIGN_WITH_BACKPOINTER
#define ALIGN_WITH_BACKPOINTER 1
#include "search_include.cpp"



template <class InputIterator>
std::string firstWord(InputIterator first, InputIterator last)
{
  return std::string(first, std::find_if(first, last, isspace));
}


// Using std::vector<char> instead of std::vector<bool> .
typedef std::vector<char> BoolVector;
BoolVector marked;
BoolVector markedGrey;
int numberOfMarked;
int remainingAtLastDbUpdate;

sequence::Fasta inputFasta;
sequence::Fasta predeterminedCentersFasta;

bool lessSpectrumIndexPair(const SpectrumIndexPair &x, const SpectrumIndexPair &y, KmerNumber dimension)
{
  return x.first[dimension] < y.first[dimension];
}

typedef std::vector<KmerCount> KmerCountVector;

template <class RandomAccessIterator>
RandomAccessIterator median(RandomAccessIterator first, RandomAccessIterator last)
{
  RandomAccessIterator mid = first + ((last - first) / 2);
  std::nth_element(first, mid, last);
  return mid;
}

typedef const KmerCount * const ConstKmerCountArray;
typedef KmerCount * const KmerCountArray;

Index twoMeans(const Index l, const Index r)
{

#if MY_DEBUG
  // debug.
  std::cerr << "\n\ttwoMeans(" << l << ",\t" << r  << ")\n" << std::flush;
#endif // MY_DEBUG

#if OLD_CODE
  std::sort(spectrumIndexes.begin() + l, spectrumIndexes.begin() + r);
#endif // OLD_CODE

#if OLD_CODE
  typedef std::vector<float> FloatVector;
  FloatVector center1(spectrumIndexes[l].first.begin(), spectrumIndexes[l].first.end());
  FloatVector center2(center1.size());
#endif // OLD_CODE

  KmerCountVector center1(spectrumIndexes[l].first.begin(), spectrumIndexes[l].first.end());
  KmerCountVector center2(center1.size());


  {
    Index i = l + 1;
    while (spectrumIndexes[l].first == spectrumIndexes[i].first)
      ++i;
    
    std::copy(spectrumIndexes[i].first.begin(), spectrumIndexes[i].first.end(), center2.begin());
  }

  assert(center1 != center2);

  typedef std::vector<int> IntVector;
  IntVector clusterNumber(r - l, 0);
  clusterNumber.front() = 1;
  clusterNumber.back() = 0;

  int clustersChanged = 100000000;


  // debug.
  int clusterChanges = 0;


  KmerCountVector counts1;
  KmerCountVector counts2;

  counts1.reserve(r - l);
  counts2.reserve(r - l);

  while((clustersChanged > 0) && ((clusterChanges <= 1) || (clustersChanged > (r - l) / 10))) {

#if MY_DEBUG
    // debug.
    std::cerr << '.' << std::flush;
    if ((l == 215699) && (r == 216534)) {
      for (Index i = l; i < r; ++i)
	std::cerr << (clusterNumber[i - l]);
      std::cerr << '\n';
    }
#endif // MY_DEBUG
    ++clusterChanges;


    clustersChanged = 0;
    
    for (Index i = l; i < r; ++i) {
      int d1 = 0;
      int d2 = 0;

      ConstKmerCountArray spectrum = spectrumIndexes[i].first.data();
      for (KmerNumber j = 0; j < number_of_k_mers; ++j) {
	d1 += abs(static_cast<int>(center1[j]) - static_cast<int>(spectrum[j]));
	d2 += abs(static_cast<int>(center2[j]) - static_cast<int>(spectrum[j]));
      }

      if (d1 != d2) {
	int newClusterNumber = (d1 < d2) ? 1 : 0;
	if (newClusterNumber != clusterNumber[i - l]) {
	  clusterNumber[i - l] = newClusterNumber;
	  ++clustersChanged;
	} // if
      }
    } // for

#if OLD_CODE

    {

      int sizeOfCluster1 = 0;
      int sizeOfCluster2 = 0;

      for (Index i = l; i < r; ++i) 
	if (clusterNumber[i - l] == 1) 
	  ++sizeOfCluster1;
	else 
	  ++sizeOfCluster2;


      if (sizeOfCluster1 == 0) {
	
	// debug.
	std::cerr << "sizeOfCluster1 == 0\n" << std::flush;

	float maxD2 = 0;
	Index maxIndex = -1;
	for (Index i = l; i < r; ++i) {
	  float d2 = 0;
	  for (KmerNumber j = 0; j < NUMBER_OF_K_MERS; ++j) {
	    d2 += fabsf(center2[j] - spectrumIndexes[i].first[j]);
	  }
	  if (d2 > maxD2) {
	    maxD2 = d2;
	    maxIndex = i;
	  }
	} // for

	assert(l <= maxIndex);
	assert(maxIndex < r);

	clusterNumber[maxIndex - l] = 1;
	clustersChanged = true;
      }
      if (sizeOfCluster2 == 0) {

	// debug.
	std::cerr << "sizeOfCluster2 == 0\n" << std::flush;

	float maxD1 = 0;
	Index maxIndex = -1;
	for (Index i = l; i < r; ++i) {
	  float d1 = 0;
	  for (KmerNumber j = 0; j < NUMBER_OF_K_MERS; ++j) {
	    d1 += fabsf(center1[j] - spectrumIndexes[i].first[j]);
	  }
	  if (d1 > maxD1) {
	    maxD1 = d1;
	    maxIndex = i;
	  }
	} // for

	assert(l <= maxIndex);
	assert(maxIndex < r);

	clusterNumber[maxIndex - l] = 0;
	clustersChanged = true;
      }
	

    }
#endif // OLD_CODE
 

    if (clustersChanged) { // Find new centers.

#if OLD_CODE
      int sizeOfCluster1 = 0;
      int sizeOfCluster2 = 0;

      std::fill(center1.begin(), center1.end(), 0);
      std::fill(center2.begin(), center2.end(), 0);

      for (Index i = l; i < r; ++i) {
	if (clusterNumber[i - l] == 1) {
	  ++sizeOfCluster1;
	  for (KmerNumber j = 0; j < NUMBER_OF_K_MERS; ++j)
	    center1[j] += spectrumIndexes[i].first[j];
	} else {
	  ++sizeOfCluster2;
	  for (KmerNumber j = 0; j < NUMBER_OF_K_MERS; ++j)
	    center2[j] += spectrumIndexes[i].first[j];
	}

      }

      assert(sizeOfCluster1 > 0);
      assert(sizeOfCluster2 > 0);

      for (KmerNumber j = 0; j < NUMBER_OF_K_MERS; ++j) {
	center1[j] /= sizeOfCluster1;
	center2[j] /= sizeOfCluster2;
      } // for
#endif // OLD_CODE


#if MY_DEBUG
      std::cerr << '(' << std::flush;
#endif // MY_DEBUG

      for (KmerNumber j = 0; j < number_of_k_mers; ++j) {

	counts1.clear();
	counts2.clear();

	for (Index i = l; i < r; ++i) {
	  if (clusterNumber[i - l] == 1) {
	    counts1.push_back(spectrumIndexes[i].first[j]);
	  } else {
	    counts2.push_back(spectrumIndexes[i].first[j]);
	  } // if (clusterNumber[i - l] == 1)
	} // for (Index i

#if MY_DEBUG
	assert(counts1.size() > 0);
	assert(counts2.size() > 0);
#endif // MY_DEBUG

	center1[j] = *median(counts1.begin(), counts1.end());
	center2[j] = *median(counts2.begin(), counts2.end());

      } // for (KmerNumber j

#if MY_DEBUG
      std::cerr << ')' << std::flush;
#endif // MY_DEBUG

    } // Find new center.


  } // while

#if MY_DEBUG
  std::cerr << 'S' << std::flush;
#endif // MY_DEBUG

  { // Split.
    Index i = l;
    Index j = r - 1;

    KmerCountVector temp(number_of_k_mers);
    
    while (j > i) {
      while (clusterNumber[i - l] == 1)
	++i;
      while (clusterNumber[j - l] == 0)
	--j;
      if (j > i) {
	std::swap(spectrumIndexes[i], spectrumIndexes[j]);
	/*
	{
	  // std::swap(spectrumIndexes[i].first, spectrumIndexes[j].first);
	  std::copy(spectrumIndexes[i].first.begin(), spectrumIndexes[i].first.end(), temp.begin());
	  std::copy(spectrumIndexes[j].first.begin(), spectrumIndexes[j].first.end(), spectrumIndexes[i].first.begin());
	  std::copy(temp.begin(), temp.end(), spectrumIndexes[i].first.begin());

	  std::swap(spectrumIndexes[i].second, spectrumIndexes[j].second);
	}
	*/
	std::swap(clusterNumber[i - l], clusterNumber[j - l]);
	++i;
	--j;
      } // if
    } // while

    // debug.
    
  } // Split.

#if MY_DEBUG
  std::cerr << 'R' << std::flush;
#endif // MY_DEBUG

  Index result = l;
  while (clusterNumber[result - l] == 1)
    ++result;

  for (Index i = result; i < r; ++i)
    assert(clusterNumber[i - l] == 0);

#if MY_DEBUG
  // debug.
  std::cerr << "clusterChanges=" << clusterChanges << ". done.\n" << std::flush;
#endif // MY_DEBUG
#if MY_DEBUG
  std::cerr << 'Q' << std::flush;
#endif // MY_DEBUG

  return result;
  
  
} // twoMeans


namespace cluster_information
{
  Index *l = 0;
  Index *r = 0;
  Index *sizeOfLeftCluster = 0;


  Index *numberOfKmers = 0;
  KmerCount **allMinCounts = 0;
  KmerCount **allMaxCounts = 0;

  int **counts = 0;
}



const KmerCount MAX_COUNT = std::numeric_limits<KmerCount>::max();
const KmerCount MIN_COUNT = 0;

const Index NULL_INDEX = -1;


void clusterSort(const Index l, const Index r, const Index index)
{

#if MY_DEBUG
  if (r - l > 10000)
    std::cerr << "clusterSort(" << l << ",\t" << r << ") " << std::flush;
#endif // MY_DEBUG

#if OLD_CODE
  ClusterInformation &info = clusters_new[index];

  clusters_new[index].l = l;
  clusters_new[index].r = r;

  KmerCountArray minCounts = clusters_new[index].minCounts;
  KmerCountArray maxCounts = clusters_new[index].maxCounts;

  std::fill_n(minCounts, number_of_k_mers, MAX_COUNT);
  std::fill_n(maxCounts, number_of_k_mers, MIN_COUNT);

  for (Index i = l; i < r; ++i)
    for (KmerNumber j = 0; j < number_of_k_mers; ++j) {
      KmerCount count = spectrumIndexes[i].first[j];

      if (count < minCounts[j])
	minCounts[j] = count;
      if (count > maxCounts[j])
	maxCounts[j] = count;
    }


#endif // OLD_CODE

#if MY_DEBUG
  std::cerr << '1' << std::flush;
#endif // MY_DEBUG

  {
    // Set 'cluster_infomation' of this cluster.
    // \for_all k \in {1...k_mer_length}: cluster_information::allMinCounts[k][index * cluster_information::numberOfKmers[k] ... index * cluster_information::numberOfKmers[k] + cluster_information::numberOfKmers[k] - 1]
    

    for (int k = 0; k < k_mer_length; ++k) {
      std::fill_n(&cluster_information::allMinCounts[k][index * cluster_information::numberOfKmers[k]], cluster_information::numberOfKmers[k], MAX_COUNT);
      std::fill_n(&cluster_information::allMaxCounts[k][index * cluster_information::numberOfKmers[k]], cluster_information::numberOfKmers[k], MIN_COUNT);
    }

#if MY_DEBUG
    std::cerr << '2' << std::flush;
#endif // MY_DEBUG

    for (Index i = l; i < r; ++i) {


      for (int k = 0; k < k_mer_length; ++k)
	std::fill_n(cluster_information::counts[k], cluster_information::numberOfKmers[k], 0);
    
      std::copy (spectrumIndexes[i].first.begin(), spectrumIndexes[i].first.end(), cluster_information::counts[k_mer_length - 1]);

      for (int k = k_mer_length - 1; k > 0; --k)
	for (int j = 0; j < cluster_information::numberOfKmers[k]; ++j) {
	  int j_prime = j / static_cast<int>(NUMBER_OF_NUCLEOTIDES);
	  cluster_information::counts[k - 1][j_prime] += cluster_information::counts[k][j];
	}

      for (int k = 0; k < k_mer_length; ++k)
	for (int j = 0; j < cluster_information::numberOfKmers[k]; ++j) {
	  KmerCount count = MAX_COUNT;
	  if (cluster_information::counts[k][j] < MAX_COUNT)
	    count = cluster_information::counts[k][j];

	  KmerCount &minCount = cluster_information::allMinCounts[k][index * cluster_information::numberOfKmers[k] + j];
	  KmerCount &maxCount = cluster_information::allMaxCounts[k][index * cluster_information::numberOfKmers[k] + j];

	  if (minCount > count)
	    minCount = count;
	  if (maxCount < count)
	    maxCount = count;
	}
    } // for (Index i = l; i < r; ++i)


    cluster_information::l[index] = l;
    cluster_information::r[index] = r;

  }

#if MY_DEBUG
  {

    std::cerr << '\n';

    for (int k = 0; k < k_mer_length; ++k) {
      for (int j = 0; j < cluster_information::numberOfKmers[k]; ++j) {
	std::cerr << static_cast<int>(cluster_information::allMinCounts[k][index * cluster_information::numberOfKmers[k] + j]) << '\t';
      }
      std::cerr << '\n';
      for (int j = 0; j < cluster_information::numberOfKmers[k]; ++j) {
	std::cerr << static_cast<int>(cluster_information::allMaxCounts[k][index * cluster_information::numberOfKmers[k] + j]) << '\t';
      }
      std::cerr << '\n';
      

    }

    std::cerr << std::flush;

  }

  std::cerr << '3' << std::flush;
#endif // MY_DEBUG


#if MEDIAN_COUNTS
  // Calculate 'medianCounts' and 'averageDeviations'.
  KmerCountArray medianCounts = info.medianCounts;
  {
    //    KmerCountArray averageDeviations = info.averageDeviations;

    KmerCountVector counts(r - l);
    KmerCount count;

    for (KmerNumber j = 0; j < NUMBER_OF_K_MERS; ++j) {
      for (Index i = l; i < r; ++i) {
	count = spectrumIndexes[i].first[j];
	counts[i - l] = count;
      } // for (Index i = l; i < r; ++i)
    
      medianCounts[j] = *median(counts.begin(), counts.end());

#if OLD_CODE
      int deviation = 0;

      for (Index i = l; i < r; ++i) {
	count = spectrumIndexes[i].first[j];
	deviation += abs(medianCounts[j] - count);
      } // for (Index i = l; i < r; ++i)

      averageDeviations[j] = deviation / (r - l);
#endif // OLD_CODE

    } // for (KmerNumber k = 0; k < NUMBER_OF_K_MERS; ++k)
  } // Calculate 'medianCounts' and 'averageDeviations'.
#endif // MEDIAN_COUNTS

#if OLD_CODE
  // Store indexes of different new counts.
  info.numberOfDifferentCounts = 0;
  for (KmerNumber j = 0; j < number_of_k_mers; ++j) 
    if ((minCounts[j] > prevMinCounts[j]) || 
	(maxCounts[j] < prevMaxCounts[j])

#if MEDIAN_COUNTS
	||
	(medianCounts[j] != prevMedianCounts[j])
#endif // MEDIAN_COUNTS

	) {
      info.indexOfDifferentCounts[info.numberOfDifferentCounts] = j;
      ++info.numberOfDifferentCounts;
    }

#endif // OLD_CODE

#if MY_DEBUG
  std::cerr << "l=" << l << "\tr=" << r << "\ti=" << index << '\n' << std::flush;
  for (KmerNumber j = 0; j < NUMBER_OF_K_MERS; ++j)
    std::cerr << static_cast<int>(minCounts[j]) << '\t';
  std::cerr << '\n' << std::flush;
  for (KmerNumber j = 0; j < NUMBER_OF_K_MERS; ++j)
    std::cerr << static_cast<int>(maxCounts[j]) << '\t';
  std::cerr << '\n' << std::flush;
#endif // MY_DEBUG

  bool allEqual = true;

  {
    for (KmerNumber j = 0; j < number_of_k_mers; ++j) {

      KmerCount &minCount = cluster_information::allMinCounts[k_mer_length - 1][index * cluster_information::numberOfKmers[k_mer_length - 1] + j];
      KmerCount &maxCount = cluster_information::allMaxCounts[k_mer_length - 1][index * cluster_information::numberOfKmers[k_mer_length - 1] + j];

      if (minCount < maxCount)
	allEqual = false;
    }
  }


  if (allEqual) {

#if OLD_CODE
    clusters_new[index].sizeOfLeftCluster = NULL_INDEX;
#endif // OLD_CODE

    cluster_information::sizeOfLeftCluster[index] = NULL_INDEX;

#if MY_DEBUG
    std::cerr << ".\n" << std::flush;
#endif // MY_DEBUG

    return;
  }

#if MY_DEBUG
  std::cerr << '4' << std::flush;
#endif // MY_DEBUG

  Index mid = twoMeans(l, r);

#if MY_DEBUG
  std::cerr << ".\n" << std::flush;
#endif // MY_DEBUG

#if OLD_CODE
  clusters_new[index].sizeOfLeftCluster = mid - l;
#endif // OLD_CODE

  cluster_information::sizeOfLeftCluster[index] = mid - l;

  clusterSort(l, mid, index + 1);
  clusterSort(mid, r, index + 2 * (mid - l));
}

IndexVector clusterSearchResults;

#if OLD_CODE
typedef const ClusterInformation * ConstClusterInfoPointer;
#endif // OLD_CODE

typedef const KmerCount * const ConstKmerCountArray;

typedef int Threshold;

Threshold maxMoreThreshold;
Threshold minMoreThreshold;



namespace layered {

  int layeredMaxMoreThresholds[MAXIMUM_K_MER_LENGTH];
  int layeredMinMoreThresholds[MAXIMUM_K_MER_LENGTH];

  int **allQueryCounts = 0;

  int layeredClusterSearch(Index index, int k)
  {
    
#if MY_DEBUG
    std::cerr << "layeredClusterSearch(" << index << ", " << k << ")\n" << std::flush;
#endif // MY_DEBUG

    int maxMore = 0;
    int minMore = 0;

    int numberOfKmers = cluster_information::numberOfKmers[k];
    KmerCount *minCounts = &cluster_information::allMinCounts[k][index * numberOfKmers];
    KmerCount *maxCounts = &cluster_information::allMaxCounts[k][index * numberOfKmers];

    int minMoreThreshold = layeredMinMoreThresholds[k];
    int maxMoreThreshold = layeredMaxMoreThresholds[k];


    int *queryCounts = allQueryCounts[k];

#if MY_DEBUG
    for (Index i = 0; i < numberOfKmers; ++i) 
      std::cerr << queryCounts[i] << '\t';
    std::cerr << '\n';
    for (Index i = 0; i < numberOfKmers; ++i) 
      std::cerr << static_cast<int>(minCounts[i]) << '\t';
    std::cerr << '\n';
    for (Index i = 0; i < numberOfKmers; ++i) 
      std::cerr << static_cast<int>(maxCounts[i]) << '\t';
    std::cerr << '\n' << std::flush;
#endif // MY_DEBUG

    {
      // Compute 'minMore' and 'maxMore'.

      int minCount;
      int maxCount;
      int queryCount;

      for (Index i = 0; i < numberOfKmers; ++i) {

	maxCount = maxCounts[i];
	queryCount = queryCounts[i];
	
	if (maxCount > queryCount) {
	  maxMore += maxCount - queryCount;
	  
	  minCount = minCounts[i];
	  if (minCount > queryCount) { 
	    minMore += minCount - queryCount;

	    if (minMore > minMoreThreshold) {
#if MY_DEBUG								
	      std::cerr << "minMore=" << minMore << " > " << layeredMinMoreThresholds[k] << "=layeredMinMoreThresholds[k]\n" << std::flush;
#endif // MY_DEBUG				
	      return 0;
	    } // if (minMore > minMoreThreshold)


	    // Second type of approximate filter
	    if (approximateFilter2 && minMore * (cluster_information::r[index] - cluster_information::l[index]) > minMoreThreshold)
	      return 0;

	  } // if (minCount > queryCount)

	} // if (maxMore > queryCount)
	  
#if OLD_CODE
	if (minCounts[i] > queryCounts[i])
	  minMore += static_cast<int>(minCounts[i]) - queryCounts[i];
	if (maxCounts[i] > queryCounts[i])
	  maxMore += static_cast<int>(maxCounts[i]) - queryCounts[i];
#endif // OLD_CODE
      } // for (Index i = 0; i < numberOfKmers; ++i)
    }

    if (maxMore <= maxMoreThreshold) {
#if MY_DEBUG
      std::cerr << "maxMore=" << maxMore << " <= " << layeredMaxMoreThresholds[k] << "=layeredMaxMoreThresholds[k]\n" << std::flush;
#endif // MY_DEBUG
      if (k + 1 < k_mer_length) {
	return layeredClusterSearch(index, k + 1);
      }

#if MY_DEBUG
      std::cerr << "Return results:\n" << std::flush;
#endif // MY_DEBUG

      for (Index i = cluster_information::l[index]; i < cluster_information::r[index]; ++i)
	clusterSearchResults.push_back(i);
      return cluster_information::r[index] - cluster_information::l[index];
    }


    return (layeredClusterSearch(index + 1, k) +
	    layeredClusterSearch(index + 2 * cluster_information::sizeOfLeftCluster[index], k));

  }

} // namespace layered


void *threaded_search_without_backpointer(void *read_args) {
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


void * threadWorker(void *context) {
  int id = (static_cast<struct Context *>(context))->id;

  pthread_mutex_lock(&searchWorkersMutexes[id]);
  while (1) {
    // Then wait for more work.
    pthread_cond_wait(&searchWorkersConds[id], &searchWorkersMutexes[id]);

    // Received work:
#if CHRIS_DEBUG
    std::cout<< "[THREAD " << id << "] received work.  Table index: " << searchWorkersArgs[id].tableIndex << "\n";
    flush(std::cout);
#endif // CHRIS_DEBUG

    
    threaded_search_without_backpointer(static_cast<void *>(&searchWorkersArgs[id]));
    pthread_mutex_lock(&allThreadsMutex);
    
    ++threads_completed;

    if (threads_completed >= running_threads) {
      pthread_cond_signal(&searchWorkersCompleted);
    }
    pthread_mutex_unlock(&allThreadsMutex);
  }

  return NULL;
}


/*
  Input: 'queryString', 'queryLength', 'queryHeader', 'queryHeaderLength' should be set.
  Database: 'sortedStrings' contains the strings to be clusterd in lexicographically sorted order. After making a cluster the sequences are removed from this vector.
  Output: Prints the cluster to STDOUT. The format of the output is determined by the value of global boolean variable 'outputMultipleAlignment'.
*/

void makeCluster()
{
#if MY_DEBUG
  std::cerr << "makeCluster()\n" << std::flush;
#endif // MY_DEBUG

  radious = static_cast<CostType>(floor((1 - similarity) * querySequenceLength));
  if (mismatches >= 0)
    radious = static_cast<CostType>(mismatches);

  if (noOverlap)
    radious *= 2;

  if (! noKmerFilter) {
    querySpectrum = countKmers(querySequencePointer, querySequenceLength);
  } // if (! noKmerFilter)



  if (! noKmerFilter) {   
    // Set allQueryCounts

    std::copy(querySpectrum.begin(), querySpectrum.end(), layered::allQueryCounts[k_mer_length - 1]);
    
    for (int k = 0; k < k_mer_length - 1; ++k)
      std::fill_n(layered::allQueryCounts[k], cluster_information::numberOfKmers[k], 0);


    for (int k = k_mer_length - 1; k > 0; --k)
      for (int j = 0; j < cluster_information::numberOfKmers[k]; ++j) {
	int j_prime = j / static_cast<int>(NUMBER_OF_NUCLEOTIDES);
	layered::allQueryCounts[k - 1][j_prime] += layered::allQueryCounts[k][j];
      }


#if MY_DEBUG
    for (int k = 0; k < k_mer_length; ++k) {
      for (int j = 0; j < cluster_information::numberOfKmers[k]; ++j)
	std::cerr << layered::allQueryCounts[k][j] << '\t';
      std::cerr << '\n' << std::flush;
    }
#endif // MY_DEBUG

    for (int k = 0; k < k_mer_length - 1; ++k)
      for (int j = 0; j < cluster_information::numberOfKmers[k]; ++j)
	if (layered::allQueryCounts[k][j] > MAX_COUNT)
	  layered::allQueryCounts[k][j] = MAX_COUNT;


  } // if (! noKmerFilter)

#if OLD_CODE
  callSplitSearch();
#endif // OLD_CODE

#if MY_DEBUG
  std::cerr << "splitSearchResults.size() == " << splitSearchResults.size() << '\n' << std::flush;
#endif // MY_DEBUG


#if MY_DEBUG
  naiveSearchResults.clear();
  
  naiveSpectrumSearch(radious * K_MER_LENGTH);

  assert(naiveSearchResults == splitSearchResults);

  std::cerr << "navieSearchResults == splitSearchResults\n" << std::flush;
#endif // MY_DEBUG

#if MY_DEBUG
  for (int i = 0; i < NUMBER_OF_K_MERS; ++i)
    std::cerr << static_cast<int>(querySpectrum[i]) << '\t';
  std::cerr << '\n' << std::flush;

  for (int j = 0; j < splitSearchResults.size(); ++j) {
    for (int i = 0; i < NUMBER_OF_K_MERS; ++i)
      std::cerr << static_cast<int>(spectrumIndexes[splitSearchResults[j]].first[i]) << '\t';
    std::cerr << '\n';
  }
#endif // MY_DEBUG

  clusterSearchResults.clear();

  if (! noKmerFilter) {
    // Set k-mer filter thresholds
    {
      for (int i = 0; i < k_mer_length; ++i)
	if (radious == 0) {
	  layered::layeredMaxMoreThresholds[i] = 0;
	  layered::layeredMinMoreThresholds[i] = 0;
	} else {
	  // radious * sqrt(i + 1)
	  // Approximation minMoreThreshold <= maxMoreThreshold
	  layered::layeredMaxMoreThresholds[i] = std::max(radious, static_cast<CostType>(querySequenceLength * 0.04)) * (i + 1);
#if OLD_CODE
	  layered::layeredMinMoreThresholds[i] = (approximateFilter ? static_cast<int>(radious * exp(approximationFactor * log(i + 1)))  : (radious * (i + 1)));
#endif // OLD_CODE
	  layered::layeredMinMoreThresholds[i] = radious * (i + 1);
	}
    } 

#if MY_DEBUG
    std::cerr << "Before layeredClusterSearch(0, 0)\n" << std::flush;
#endif //MY_DEBUG

    layered::layeredClusterSearch(0, 0);

#if MY_DEBUG
    for (int i = 0; i < static_cast<int>(clusterSearchResults.size()); ++i)
      std::cerr << clusterSearchResults[i] << '\t';
    std::cerr << '\n' << std::flush;

    exit(EXIT_FAILURE);
#endif // MY_DEBUG


#if NEW_CODE
    clusterSearch(clusters.data());
#endif // NEW_CODE
  } // if (! noKmerFilter)

  IndexVector unmarkeds;


  if (! noKmerFilter) { 
    // Remove already marked sequences.
    unmarkeds.reserve(clusterSearchResults.size());

    for(IndexVector::const_iterator it = clusterSearchResults.begin(); it != clusterSearchResults.end(); ++it){
      Index lexicographicIndex = spectrumIndexes[*it].second;
      if (!marked[lexicographicIndex])
	unmarkeds.push_back(lexicographicIndex);
    }
    const Index numberOfUnmarkeds = unmarkeds.size();
  
    std::sort(unmarkeds.begin(), unmarkeds.end());

#if MY_DEBUG							
    std::cout << "Query:\t" << querySequencePointer << '\n';
  
    for (Index i = 0; i < numberOfUnmarkeds; ++i)
      std::cout << lexicographicSequencePointers[unmarkeds[i]] << '\n';
#endif // MY_DEBUG

    sortedStrings.resize(numberOfUnmarkeds);

    PositionVector sortedStringsLengths(sortedStrings.size());
  
    for (Index i = 0; i < numberOfUnmarkeds; ++i) {
      sortedStrings[i] = lexicographicSequencePointers[unmarkeds[i]];
      sortedStringsLengths[i] = sequenceLengths[unmarkeds[i]];
    }

    searchResults.clear();

    for(uint16_t i = 0; i < num_threads; ++i){
      word[i][1] = invalidChar;
      wordMinLength[i][1] = MAX_LENGTH + 1;
    }

    
    size_t num_threads_to_use = num_threads;
    int interval_size = static_cast<int>(ceil(sortedStrings.size() / static_cast<double>(num_threads)));//sortedStrings.size() / num_threads_to_use;

    if (sortedStrings.size() < num_threads) {
      num_threads_to_use = 1;
      interval_size = sortedStrings.size();
    }

#if CHRIS_DEBUG
      std::cout << "Size: " <<  sortedStrings.size() << "\n";
      //std::cout << "Building first trees for table: " << i << " from: " <<  << " to " << end << ", interval_size " << chunkIndexes[i].second << "\n" << std::flush;
#endif // CHRIS_DEBUG

    running_threads = 0;
    for (unsigned int i = 0, currentTable = 0; i < sortedStrings.size(); i += interval_size, ++currentTable) {
      unsigned int end = i + interval_size;
      if (i + interval_size > sortedStrings.size()) {
        end = sortedStrings.size();
      }
      else {
        end = i + interval_size;
      }

      sortedStringsLengthsMinTree[currentTable].resize(2 * (end - i));
      sortedStringsLengthsMaxTree[currentTable].resize(2 * (end - i));

#if CHRIS_DEBUG
        std::cout << "Kmer: Building second trees from: " << i << " to " << end << "\n" << std::flush; 
        std::cout << "Kmer: Size of sorted string min/max: " << sortedStringsLengthsMinTree[currentTable].size() << ", " 
            << sortedStringsLengthsMaxTree[currentTable].size() << "\n" << std::flush;
#endif // CHRIS_DEBUG

      buildMinTree(sortedStringsLengths.begin() + i, sortedStringsLengths.begin() + end, sortedStringsLengthsMinTree[currentTable].begin());
      buildMaxTree(sortedStringsLengths.begin() + i, sortedStringsLengths.begin() + end, sortedStringsLengthsMaxTree[currentTable].begin());
    }

    if (numberOfUnmarkeds) {
      if (outputMultipleAlignment)
	search(0, numberOfUnmarkeds - 1, 0, 0, 0, 0);
      else {
      
      running_threads = 0;
      threads_completed = 0;

#if CHRIS_DEBUG
        std::cout << "Size: " << sortedStrings.size() << ", numberOfUnmarkeds: " << numberOfUnmarkeds 
        << ", interval_size" << interval_size<< "\n" <<std::flush;
#endif // CHRIS_DEBUG

      for (int i = 0; i < numberOfUnmarkeds; i += interval_size) {
        ++running_threads;
        pthread_mutex_lock(&searchWorkersMutexes[running_threads-1]);

        struct thread_args *args = &searchWorkersArgs[running_threads-1];
        args->l = 0;

        if (i + interval_size > numberOfUnmarkeds) {
          args->r = (numberOfUnmarkeds - 1) - i;
        }
        else {
          args->r = interval_size - 1;
        }
        args->pos = 0;
        args->intervalIndex = 0;
        args->tableIndex = running_threads-1;
        args->globalOffset = i;

#if CHRIS_DEBUG
        std::cout << "  Assigning work to thread: " << i << ", l: " << args->l
            << ", r: " << args->r << "\n" <<std::flush;
#endif // CHRIS_DEBUG

        pthread_cond_signal(&searchWorkersConds[running_threads-1]);
        pthread_mutex_unlock(&searchWorkersMutexes[running_threads-1]);

      }

      // Wait for all threads to complete.
      pthread_mutex_lock(&allThreadsMutex);
      while (threads_completed < running_threads) {
        pthread_cond_wait(&searchWorkersCompleted, &allThreadsMutex);
      }
      pthread_mutex_unlock(&allThreadsMutex);

      for(int i = 0; i < running_threads; i++)  {
        searchResults.insert(searchResults.end(), threadSearchResults[i].begin(), threadSearchResults[i].end());
        threadSearchResults[i].clear();
      }
    }

	//search_without_backpointer(0, numberOfUnmarkeds - 1, 0, 0, 0, 0);
    }

  } // if (! noKmerFilter)
  else {
    // No k-mer filter

    for(uint16_t i = 0; i < num_threads; ++i){
      word[i][1] = invalidChar;
      wordMinLength[i][1] = MAX_LENGTH + 1;
    }
    searchResults.clear();

    if (outputMultipleAlignment)
      search(0, sortedStrings.size() - 1, 0, 0, 0, 0);
    else {
      int num_threads_to_use = num_threads;
      int interval_size = static_cast<int>(ceil(sortedStrings.size() / static_cast<double>(num_threads)));
      
      if (sortedStrings.size() < num_threads) {
        num_threads_to_use = 1;
        interval_size = sortedStrings.size();
      }

      running_threads = 0;
      threads_completed = 0;

      for (unsigned int i = 0; i < sortedStrings.size(); i += interval_size) {
        ++running_threads;

        // Create the thread variables.
        pthread_mutex_lock(&searchWorkersMutexes[running_threads - 1]);

        struct thread_args *args = &searchWorkersArgs[running_threads - 1];
        args->l = 0;

        if (i + interval_size > sortedStrings.size()) {
          args->r = (sortedStrings.size() - 1) - i;
        }
        else {
          args->r = interval_size - 1;
        }
        args->pos = 0;
        args->intervalIndex = 0;
        args->tableIndex = running_threads-1;
        args->globalOffset = i;

#if CHRIS_DEBUG
        std::cout << "Assigning work to thread: " << i << ", l: " << args->l
            << ", r: " << args->r << "\n" <<std::flush;
#endif // CHRIS_DEBUG


        pthread_cond_signal(&searchWorkersConds[running_threads-1]);
        pthread_mutex_unlock(&searchWorkersMutexes[running_threads-1]);
      }

      // Wait for all threads to complete.
      pthread_mutex_lock(&allThreadsMutex);
      while (threads_completed < running_threads) {
        pthread_cond_wait(&searchWorkersCompleted, &allThreadsMutex);
      }
      pthread_mutex_unlock(&allThreadsMutex);

      for(int i = 0; i < running_threads; i++)  {
        searchResults.insert(searchResults.end(), threadSearchResults[i].begin(), threadSearchResults[i].end());
        threadSearchResults[i].clear();
      }

      //search_without_backpointer(0, sortedStrings.size() - 1, 0, 0);
    }
  }


  { // Output cluster alignment.
	

    // Output cluster size.
    if (outputMultipleAlignment) {
      std::cout << CLUSTER_INDICATOR << searchResults.size() + 1 << '\n';
    }


    typedef std::vector<PositionType> CountsVector;


    // This vector will contain the maximum number of gaps that happen after position [i] in query string in ANY of the alignments to search results.
    CountsVector gapCounts(querySequenceLength);
	
	
    if (outputMultipleAlignment) { // Compute gapCounts[0..].

      for (SearchResultVector::const_iterator resultIterator = searchResults.begin(); resultIterator != searchResults.end(); ++resultIterator) {
	for (PositionType i = 0; i < static_cast<PositionType>(resultIterator->func.size() - 1); ++i)
	  gapCounts[i] = std::max(gapCounts[i], resultIterator->func[i + 1] - resultIterator->func[i] - 1);

      } // for
    } // if (outputMultipleAlignment)

    // This string will contain the query sequence with appropriate number of gaps inserted at each position to produce a multiple alignment of the cluster.
    std::string clusterCenterSequenceWithGaps;

    if (outputMultipleAlignment) { // Compute clusterCenterSequenceWithGaps.
      for (PositionType positionInQuery = 0; positionInQuery < querySequenceLength; ++positionInQuery) {
	clusterCenterSequenceWithGaps.append(gapCounts[positionInQuery], gapChar);
	clusterCenterSequenceWithGaps.push_back(querySequencePointer[positionInQuery]);
      } // for
    } // if (outputMultipleAlignment)

    // Output cluster center sequence with gaps from multiple alignment.
    if (outputMultipleAlignment) {
      std::cout << '>' << queryHeaderPointer << '\n'
		<< clusterCenterSequenceWithGaps << '\n';
    } else {
      if (useFullQueryHeader) {
        if (!assignAmbiguous) {
          std::cout << queryHeaderPointer << '\t';
        } else {
          if (!printInvertedIndex)
            center_nonambig_counts[queryHeaderPointer] = 1;
        }
      } else {
        if (!assignAmbiguous) {
          std::cout << firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength) << '\t';
        } else {
          if (!printInvertedIndex)
            center_nonambig_counts[firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength)] = 1;
        }
      }
    } 
	
    // Output multiple alignment of sequences in current cluster.
    for (SearchResultVector::const_iterator resultIterator = searchResults.begin(); resultIterator != searchResults.end(); ++resultIterator) {

      Index lexicographicIndex;
      if (! noKmerFilter) {
	lexicographicIndex = unmarkeds[resultIterator->number];
      } else {
	lexicographicIndex = lexicographicIndexOfSortedStringsIndexes[resultIterator->number];
      }
      ConstCharPointer s = sortedStrings[resultIterator->number];

      if (!assignAmbiguous)
        markedGrey[lexicographicIndex] = true;

      if ((!noOverlap) || resultIterator->cost < radious / 2) {
      
	if (!marked[lexicographicIndex]) {
    if (!assignAmbiguous) {
  	  marked[lexicographicIndex] = true;
  	  ++numberOfMarked;
    }

	  ConstCharPointer headerPointer = headerPointers[lexicographicIndex];
	  PositionType headerLength = headerLengths[lexicographicIndex];

	  if (outputMultipleAlignment) {
	    std::string alignedSequence; 
	    for (PositionType i = 0; i < static_cast<PositionType>(gapCounts.size()); ++i) {
	      for (PositionType j = resultIterator->func[i + 1] - resultIterator->func[i]; j < gapCounts[i] + 1; ++j)
		alignedSequence.push_back(gapChar);
	      for (PositionType j = resultIterator->func[i]; j < resultIterator->func[i + 1]; ++j)
		alignedSequence.push_back(s[j]);
	    } // for (PositionType i = 0; i < gapCounts.size() - 1; ++i)

	    std::cout << '>' << headerPointer << '\n'
		      << alignedSequence << '\n';
	  } else {
      if (useFullQueryHeader) {
        if (!assignAmbiguous) {
          std::cout << headerPointer << '\t';
        }
        else {
          if (printInvertedIndex) {
            std::cout << headerPointer << '\t' << queryHeaderPointer << '\n';
          } else {
            inverted_index[headerPointer].push_back(queryHeaderPointer);
            if (inverted_index[headerPointer].size() == 1) {
              center_nonambig_counts[queryHeaderPointer] += 1;
            } else if (inverted_index[headerPointer].size() == 2) {
              center_nonambig_counts[inverted_index[headerPointer][0]] -= 1;
            }
          }
        }

      } else {
        if (!assignAmbiguous) {
          std::cout << firstWord(headerPointer, headerPointer + headerLength) << '\t';
        } else {
          if (printInvertedIndex) {
            std::cout << firstWord(headerPointer, headerPointer + headerLength) << '\t' 
                << firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength) << '\n';
          } else {
            inverted_index[firstWord(headerPointer, headerPointer + headerLength)].push_back(
                firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength));

            //inverted_index[firstWord(headerPointer, headerPointer + headerLength)].push_back(
            //  firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength));
            if (inverted_index[firstWord(headerPointer, headerPointer + headerLength)].size() == 1) {
              center_nonambig_counts[firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength)] += 1;
            } else if (inverted_index[firstWord(headerPointer, headerPointer + headerLength)].size() == 2) {
              //std::cout << inverted_index[firstWord(headerPointer, headerPointer + headerLength)];
              //std::cout << "Prev index: " << inverted_index[firstWord(headerPointer, headerPointer + headerLength)][0] << "+";
              center_nonambig_counts[inverted_index[firstWord(headerPointer, headerPointer + headerLength)][0]] -= 1;
            }
          }
        }
      }

	  }
	    
	} // if (!marked[lexicographicIndex])

      } // if ((!noOverlap) || resultIterator->cost < radious / 2)

    } // for
	
    if (!assignAmbiguous)
      std::cout << '\n';

  } // Output cluster alignment.

}


int main(int argc, char *argv[])
{

  { // Parse command line options.
    
    const std::string PROGRAM_NAME = "dnaclust";
    const std::string USAGE_DESCRIPTION =
      "  The output is written to STDOUT."
      " Each line will contain the ids of the sequences in each cluster, and the first id of each line is the cluster representative.\n"
      "  Example: To cluster a set of 16S rRNA fragments at 0.98 similarity use:\n"
      "  ./dnaclust file.fasta -l -s 0.98 > clusters \n"
      "  You can optionally specify a k-mer length for the filter."
      " The longer k-mers use more memory. "
      " Also the filter will be more specific with longer k-mers."
      " The default log_4(median length) should be good for most cases.\n"
      ;
#if OLD_CODE
	  "  For the fastest running time please specify k-mer length around log_4(sequence-length).\n"
	  "  Example: To cluster a set of 16S RNA fragments of length around 250 at 0.98 similarity use.\n"
	  "  ./dnaclust file.fasta -l -s 0.98 -k 3 > clusters \n"
#endif // OLD_CODE


    po::options_description desc("Options");

    std::string inputFileName;
    std::string predeterminedCentersFileName;
    
    desc.add_options()
      ("help,h", "produce help message")
      ("similarity,s", po::value<float>(&similarity)->default_value(DEFAULT_SIMILARITY), "set similarity between cluster center and cluster sequences")
      ("input-file,i", po::value<std::string>(&inputFileName), "input file")
      ("predetermined-cluster-centers,p", po::value<std::string>(&predeterminedCentersFileName), "file containing predetermined cluster centers")
      ("recruit-only,r", po::bool_switch(&recruit_only), "when used with 'predetermined-cluster-centers' option, only clusters the input sequences that are similar to the predetermined centers")
      //("multiple-alignment,m", po::bool_switch(&outputMultipleAlignment)->default_value(false), "produce multiple alignment for each cluster")
      ("header,d", "output header line indicating run options")
      ("left-gaps-allowed,l", "allow for gaps on the left of shorter string in semi-global alignment")
      ("k-mer-length,k", po::value<int>(&k_mer_length)/*->default_value(DEFAULT_K_MER_LENGTH) */, "length of k-mer for filtering")
#if OLD_CODE
      ("approximate-filter,a", po::bool_switch(&approximateFilter)->default_value(false), "use faster approximate k-mer filter")
      ("approximation-factor,f", po::value<double>(&approximationFactor)->default_value(0.5), "determine amount of approximation, by a floating point value between zero and one. smaller values result in more approximation.")
      ("approximate-filter2", po::bool_switch(&approximateFilter2)->default_value(false), "use second type of approximate k-mer filter")
#endif // OLC_CODE
      ("approximate-filter", po::bool_switch(&approximateFilter2)->default_value(false), "use faster approximate k-mer filter")
      ("k-mer-filter", po::bool_switch(&kmerFilter)->default_value(false), "use k-mer filter")
      ("no-k-mer-filter", po::bool_switch(&noKmerFilter)->default_value(true), "do not use k-mer filter")
      ("no-overlap", po::bool_switch(&noOverlap)->default_value(false), "cluster some of sequences such that the cluster centers are at distance at least two times the radius of the clusters")
      ("threads,t", po::value<size_t>(&num_threads), "number of threads")
      ("use-full-query-header,u", po::bool_switch(&useFullQueryHeader)->default_value(false), "use the full query header instead of the first word")
      ("mismatches,m", po::value<float>(&mismatches)->default_value(DEFAULT_MISMATCHES), "number of mismatches allowed from cluster center")
      ("assign-ambiguous,a", po::bool_switch(&assignAmbiguous), "assign ambiguous reads to clusters based on abundances of non-ambiguous reads")
      ("random-seed,e", po::value<int>(&random_seed), "Seed for random number generator")
      ("print-inverted-index", po::bool_switch(&printInvertedIndex), "Print mapping from sequence to each center")
      ;

    po::positional_options_description pod;
    pod.add("input-file", -1);

    typedef std::vector<std::string> StringVector;
    StringVector options(argv, argv + argc);
    
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
    notify(vm);
    
    if(vm.count("help") || 
       (inputFileName.empty()))
      {
	std::cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] FastaFile\n"
		  << USAGE_DESCRIPTION
		  << desc << '\n';
	exit(EXIT_FAILURE);
      }

    if(kmerFilter)
      noKmerFilter = false;

    std::ifstream inputFile(inputFileName.c_str());
    inputFile >> inputFasta;

    numberOfInputSequences = inputFasta.size();

    // Initialize all threads.
    for (size_t i = 0; i < MAX_THREADS; ++i) {
      pthread_mutex_init(&searchWorkersMutexes[i], NULL);
      pthread_cond_init(&searchWorkersConds[i], NULL);

      //struct thread_args *args = &args_array[running_threads-1];
      searchWorkersContext[i].id = i;
      pthread_create(&searchWorkers[i], NULL, threadWorker, static_cast<void *>(&searchWorkersContext[i]));
    }
    pthread_mutex_init(&allThreadsMutex, NULL);
    pthread_cond_init(&searchWorkersCompleted, NULL);
    pthread_mutex_init(&running_mutex, NULL);

    if (vm.count("predetermined-cluster-centers")) {
      std::ifstream predeterminedCentersFile(predeterminedCentersFileName.c_str());
      predeterminedCentersFile >> predeterminedCentersFasta;
    }


    if (vm.count("header")){ // Output command line options at the first line of output.
      std::cout << HEADER_INDICATOR;
      for (StringVector::const_iterator i = options.begin(); i != options.end(); ++i)
	std::cout << ' ' << *i;
      std::cout << '\n';
    }

    if (vm.count("left-gaps-allowed")) {
      leftGapsAllowed = true;
    }

    if (vm.count("k-mer-length") == 0) {
      // Try to guess appropriate k-mer length based on median sequence length.

      typedef std::vector<size_t> size_vector;
      size_vector sizes;
      for (sequence::Fasta::const_iterator i = inputFasta.begin(); i != inputFasta.end(); ++i)
	sizes.push_back(i->sequence.length());
      size_t median_size = *median(sizes.begin(), sizes.end());
      k_mer_length = static_cast<int>(floor(log(median_size) / log(NUMBER_OF_NUCLEOTIDES)));
    }

    // Initializing this at the beginning is crucial.
    number_of_k_mers = power(static_cast<int>(NUMBER_OF_NUCLEOTIDES), k_mer_length);

  }


  {

    { // Calculate global vectors indexed by 'lexicographical index': lexicographicalSequencePointers, sequenceLengths, headerPointers, headerLengths

      typedef std::pair<ConstCharPointer, Index> PointerIndexPair;
      typedef std::vector<PointerIndexPair> PointerIndexVector;
    
      PointerIndexVector pointerFastaIndexes(numberOfInputSequences);
      lexicographicSequencePointers.resize(numberOfInputSequences);

      for (Index i = 0; i < numberOfInputSequences; ++i) {
	lexicographicSequencePointers[i] = inputFasta[i].sequence.c_str();

	pointerFastaIndexes[i].first = inputFasta[i].sequence.c_str();
	pointerFastaIndexes[i].second = i;
      }

      ternary_sort::ternarySort(lexicographicSequencePointers.data(), numberOfInputSequences);

      std::sort(pointerFastaIndexes.begin(), pointerFastaIndexes.end());
    
      PointerIndexVector pointerLexicographicIndexes(numberOfInputSequences);

      for (Index i = 0; i < numberOfInputSequences; ++i) {
	pointerLexicographicIndexes[i].first = lexicographicSequencePointers[i];
	pointerLexicographicIndexes[i].second = i;
      }

      std::sort(pointerLexicographicIndexes.begin(), pointerLexicographicIndexes.end());

      IndexVector fastaIndexOfLexicographicIndexes(numberOfInputSequences);
      IndexVector lexicographicIndexOfFastaIndexes(numberOfInputSequences);

      for (Index i = 0; i < numberOfInputSequences; ++i) {
	assert(pointerFastaIndexes[i].first == pointerLexicographicIndexes[i].first);
      
	fastaIndexOfLexicographicIndexes[pointerLexicographicIndexes[i].second] = pointerFastaIndexes[i].second;
	lexicographicIndexOfFastaIndexes[pointerFastaIndexes[i].second] = pointerLexicographicIndexes[i].second;

      }

      headerPointers.resize(numberOfInputSequences);
      headerLengths.resize(numberOfInputSequences);
      sequenceLengths.resize(numberOfInputSequences);

      for (Index i = 0; i < numberOfInputSequences; ++i) {
	Index fastaIndex = fastaIndexOfLexicographicIndexes[i];
	headerPointers[i] = inputFasta[fastaIndex].header.c_str();
	headerLengths[i] = inputFasta[fastaIndex].header.length();
	sequenceLengths[i] = inputFasta[fastaIndex].sequence.length();
      }
    } // Calculate global vectors indexed by 'lexicographical index': lexicographicalSequencePointers, sequenceLengths, headerPointers, headerLengths



    typedef std::pair<PositionType, Index> LengthIndexPair;
    typedef std::vector<LengthIndexPair> LengthIndexVector;
#if MY_DEBUG
    typedef std::pair<std::pair<PositionType, ConstCharPointer>, Index> LengthIndexDebugPair;
    typedef std::vector<LengthIndexDebugPair> LengthIndexVector;
#endif // MY_DEBUG
    LengthIndexVector lengthLexicographicIndexes(numberOfInputSequences);
    
    for (Index i = 0; i < numberOfInputSequences; ++i) {

      lengthLexicographicIndexes[i].first = sequenceLengths[i];
      lengthLexicographicIndexes[i].second = i;
#if MY_DEBUG
      // Assume i is fastaIndex;
      Index lexicographicIndex = lexicographicIndexOfFastaIndexes[i];
      lengthLexicographicIndexes[i].first.first = sequenceLengths[lexicographicIndex];
      lengthLexicographicIndexes[i].first.second = lexicographicSequencePointers[lexicographicIndex];
      lengthLexicographicIndexes[i].second = lexicographicIndex;
#endif // MY_DEBUG
    }


    std::sort(lengthLexicographicIndexes.begin(), lengthLexicographicIndexes.end(), std::greater<LengthIndexPair>());
#if MY_DEBUG
    std::sort(lengthLexicographicIndexes.begin(), lengthLexicographicIndexes.end(), std::greater<LengthIndexDebugPair>());
#endif // MY_DEBUG

    if (! noKmerFilter) { 

      spectrumIndexes.resize(numberOfInputSequences);

      for (Index i = 0; i < numberOfInputSequences; ++i) {
	spectrumIndexes[i].first = countKmers(lexicographicSequencePointers[i], sequenceLengths[i]);
	spectrumIndexes[i].second = i;
      }


      
      { // Sorting the k-mer spectrums in a clever way.
	
	{ // Allocate cluster_information
	  cluster_information::l = new Index[2 * numberOfInputSequences];
	  cluster_information::r = new Index[2 * numberOfInputSequences];
	  cluster_information::sizeOfLeftCluster = new Index[2 * numberOfInputSequences];

	  cluster_information::numberOfKmers = new Index[k_mer_length];
	  for (int i = 0; i < k_mer_length; ++i)
	    cluster_information::numberOfKmers[i] = power(static_cast<int>(NUMBER_OF_NUCLEOTIDES), i + 1);

	  cluster_information::allMinCounts = new KmerCount *[k_mer_length];
	  for (int i = 0; i < k_mer_length; ++i)
	    cluster_information::allMinCounts[i] = new KmerCount[cluster_information::numberOfKmers[i] * numberOfInputSequences * 2];

	  cluster_information::allMaxCounts = new KmerCount *[k_mer_length];
	  for (int i = 0; i < k_mer_length; ++i)
	    cluster_information::allMaxCounts[i] = new KmerCount[cluster_information::numberOfKmers[i] * numberOfInputSequences * 2];

	  cluster_information::counts = new int *[k_mer_length];
	  for(int i = 0; i < k_mer_length; ++i)
	    cluster_information::counts[i] = new int[cluster_information::numberOfKmers[i]];

	}

	// Allocate allQueryCounts
	{
	  layered::allQueryCounts = new int *[k_mer_length];
	  for (int i = 0; i < k_mer_length; ++i)
	    layered::allQueryCounts[i] = new int[cluster_information::numberOfKmers[i]];
	}

#if MY_DEBUG
	std::cerr << "Before clusterSort()\n" << std::flush;
#endif // MY_DEBUG

	{
	  clusterSort(0, numberOfInputSequences, 0);


#if MY_DEBUG
	  std::cerr << "After clusterSort.\n" << std::flush;
#endif // MY_DEBUG

	  remainingAtLastDbUpdate = numberOfInputSequences;
	}

      }



#if MY_DEBUG
      for (Index i = 0; i < numberOfInputSequences; ++i) {
	std::cout << i << ":\t";
	for (KmerNumber j = 0; j < NUMBER_OF_K_MERS; ++j)
	  std::cout << static_cast<int>(spectrumIndexes[i].first[j]) << '\t';
	std::cout << '\n';
      }

    exit(EXIT_FAILURE);
#endif // MY_DEBUG

    } // if (! noKmerFilter)

    if (noKmerFilter) {
      sortedStrings = lexicographicSequencePointers;
      lexicographicIndexOfSortedStringsIndexes.resize(sortedStrings.size());
      for (Index i = 0; i < static_cast<Index>(sortedStrings.size()); ++i)
	lexicographicIndexOfSortedStringsIndexes[i] = i;

      // Go through each thread and resize all the threads
      int num_threads_to_use = num_threads;
      int interval_size = static_cast<int>(ceil(sortedStrings.size() / static_cast<double>(num_threads)));//sortedStrings.size() / num_threads_to_use;

      if (sortedStrings.size() < num_threads) {
        num_threads_to_use = 1; //sortedStrings.size(); //1; //sortedStrings.size();
        interval_size = sortedStrings.size();//1;//sortedStrings.size();
      }

      //std::cout << "Using threads: " <<  num_threads_to_use << "\n";
      //flush(std::cout);
      
      
      running_threads = 0;
      //std::cout << "Sorted strings size: " << sortedStrings.size() << "\n"; 
      //flush(std::cout);

#if CHRIS_DEBUG
       std::cout << "Size: " <<  sortedStrings.size() << "\n";
      //std::cout << "Building first trees for table: " << i << " from: " << chunkIndexes[i].first << " to " << end << ", interval_size " << chunkIndexes[i].second << "\n" << std::flush;
#endif // CHRIS_DEBUG

      for (unsigned int i = 0, currentTable = 0; i < sortedStrings.size(); i += interval_size, ++currentTable) {
        unsigned int end = i + interval_size;
        if (i + interval_size > sortedStrings.size()) {
          end = sortedStrings.size();
        }
        else {
          end = i + interval_size;
        }

        sortedStringsLengthsMinTree[currentTable].resize(2 * (end - i));
        sortedStringsLengthsMaxTree[currentTable].resize(2 * (end - i));
  
#if CHRIS_DEBUG
        std::cout << "Building first trees from: " << i << " to " << end << "\n"; 
#endif // CHRIS_DEBUG

        
        buildMinTree(sequenceLengths.begin() + i, sequenceLengths.begin() + end, sortedStringsLengthsMinTree[currentTable].begin());
        buildMaxTree(sequenceLengths.begin() + i, sequenceLengths.begin() + end, sortedStringsLengthsMaxTree[currentTable].begin());
      }

      remainingAtLastDbUpdate = numberOfInputSequences;
    }

    marked.resize(numberOfInputSequences);
    std::fill(marked.begin(), marked.end(), false);

    markedGrey.resize(numberOfInputSequences);
    std::fill(markedGrey.begin(), markedGrey.end(), false);
    
    numberOfMarked = 0;
    

    // This MUST be called before any call to any of the search functions.
    initialize();
    

    // Loop over predetermined cluster centers.
    for (sequence::Fasta::const_iterator predeterminedIterator = predeterminedCentersFasta.begin(); predeterminedIterator != predeterminedCentersFasta.end();  ++predeterminedIterator) {
      querySequencePointer = predeterminedIterator->sequence.c_str();
      querySequenceLength = predeterminedIterator->sequence.length();
      queryHeaderPointer = predeterminedIterator->header.c_str();
      queryHeaderLength = predeterminedIterator->header.length();

      makeCluster();

    } // for


    if (recruit_only) {
      // If we are assigning ambiguous reads based on abundance, handle it here.
      if (assignAmbiguous) {
        /* initialize random seed: */
        if (random_seed == 0) {
          srand(time(NULL));
        } else {
          srand(random_seed);
        }
//        std::cout << "********" << "\n";
      
        typedef std::map< std::string, std::vector<std::string> >::const_iterator MapIterator;
        for (MapIterator iter = inverted_index.begin(); iter != inverted_index.end(); iter++)
        {
            //std::cout << "Key: " << iter->first << std::endl << "Values:" << std::endl;

            if (iter->second.size() == 1) {
              final_clusters[iter->second[0]].push_back(iter->first);
            } 
            else{
              std::vector<std::string> bins;
              typedef std::vector<std::string>::const_iterator ListIterator;
              for (ListIterator center_iter = iter->second.begin(); center_iter != iter->second.end(); center_iter++) {
                  for (int z = 0; z < center_nonambig_counts[*center_iter]; ++z) {
                    bins.push_back(*center_iter);
                  }
              }

              int index = rand() % bins.size();

#if ABUN_DEBUG
              std::cout << "Bins: " << bins << ", selected index: " << index << 
                  "(" << bins[index] << ")" << "\n";
#endif // ABUN_DEBUG

              final_clusters[bins[index]].push_back(iter->first);
  //            std::cout << " " << *list_iter << std::endl;
            }
        }

        for (MapIterator iter = final_clusters.begin(); iter != final_clusters.end(); iter++) {
          std::cout << iter->first << "\t";

          typedef std::vector<std::string>::const_iterator ListIterator;
          for (ListIterator seq_iter = iter->second.begin(); seq_iter != iter->second.end(); seq_iter++) {
            std::cout << *seq_iter << "\t";
          }

          std::cout << "\n";
        }
      }

      return EXIT_SUCCESS;
    }

#if TIMING
    boost::posix_time::ptime beginning = boost::posix_time::second_clock::local_time();
#endif // TIMING

    for (Index i = 0; i < numberOfInputSequences; ++i) {

      Index lexicographicIndex = lengthLexicographicIndexes[i].second;
      if (!marked[lexicographicIndex]) {
	if (!markedGrey[lexicographicIndex]) {
	  markedGrey[lexicographicIndex] = true;
	  marked[lexicographicIndex] = true;
	  ++numberOfMarked;


	if ((numberOfInputSequences - numberOfMarked < remainingAtLastDbUpdate * 0.5) && (remainingAtLastDbUpdate > numberOfInputSequences * 0.05) && (numberOfInputSequences - numberOfMarked > 1)) {
	  if (! noKmerFilter) {

	    Index i = 0;
	    Index j = remainingAtLastDbUpdate - 1;
	  
	    while (i < j) {
	      while (!marked[spectrumIndexes[i].second])
		++i;
	      while (marked[spectrumIndexes[j].second])
		--j;
	      
	      if (i < j) {
		std::swap(spectrumIndexes[i], spectrumIndexes[j]);

		++i;
		--j;
	      } // if (i < j)
	    } // while (i < j)

	    remainingAtLastDbUpdate = numberOfInputSequences - numberOfMarked;

#if TIMING								
	    std::cerr << (boost::posix_time::second_clock::local_time() - beginning)  << '\t' << "Database update.\n" << std::flush;
#endif // TIMING
							
	    clusterSort(0, remainingAtLastDbUpdate, 0);
	    
#if TIMING								
	    std::cerr << (boost::posix_time::second_clock::local_time() - beginning)  << '\t' << "...done.\n" << std::flush;
#endif // TIMING				
	  } // if (! noKmerFilter)
	  else { // No k-mer filter

	    ConstCharPointerVector newSortedStrings;
	    PositionVector newSortedStringsLengths;
	    IndexVector newLexicographicIndexOfSortedStringsIndexes;
	  

	    for (Index i = 0; i < static_cast<Index>(sortedStrings.size()); ++i) {
	      Index lexicographicIndex = lexicographicIndexOfSortedStringsIndexes[i];
	      if (! marked[lexicographicIndex]){
		newSortedStrings.push_back(lexicographicSequencePointers[lexicographicIndex]);
		newLexicographicIndexOfSortedStringsIndexes.push_back(lexicographicIndex);
		newSortedStringsLengths.push_back(sequenceLengths[lexicographicIndex]);
	      }
	    }
					 
	    sortedStrings.swap(newSortedStrings);
	    lexicographicIndexOfSortedStringsIndexes.swap(newLexicographicIndexOfSortedStringsIndexes);
	  
      size_t num_threads_to_use = num_threads;
      int interval_size = static_cast<int>(ceil(sortedStrings.size() / static_cast<double>(num_threads)));//sortedStrings.size() / num_threads_to_use;

      if (sortedStrings.size() < num_threads) {
        num_threads_to_use = 1;//sortedStrings.size(); //1;
        interval_size = sortedStrings.size();//1; //sortedStrings.size();
      }
      
#if CHRIS_DEBUG
       std::cout << "Size: " <<  sortedStrings.size() << "\n";
      //std::cout << "Building first trees for table: " << i << " from: " << chunkIndexes[i].first << " to " << end << ", interval_size " << chunkIndexes[i].second << "\n" << std::flush;
#endif // CHRIS_DEBUG
      
      for (unsigned int i = 0, currentTable = 0; i < sortedStrings.size() ; i += interval_size, ++currentTable) {

        unsigned int end = i + interval_size;
        if (i + interval_size > sortedStrings.size()) {
          end = sortedStrings.size();
        }
        else {
          end = i + interval_size;
        }
        sortedStringsLengthsMinTree[currentTable].resize(2 * (end - i));
        sortedStringsLengthsMaxTree[currentTable].resize(2 * (end - i));

#if CHRIS_DEBUG
        std::cout << "Building second trees from: " << i << " to " << end << "\n" << std::flush; 
        std::cout << "Size of sorted string min/max: " << sortedStringsLengthsMinTree[currentTable].size() << ", " 
            << sortedStringsLengthsMaxTree[currentTable].size() << "\n" << std::flush;
        //std::cout << "Building first trees for table: " << i << " from: " << chunkIndexes[i].first << " to " << end << ", interval_size " << chunkIndexes[i].second << "\n" << std::flush;
#endif // CHRIS_DEBUG
        buildMinTree(newSortedStringsLengths.begin() + i, newSortedStringsLengths.begin() + end, sortedStringsLengthsMinTree[currentTable].begin());
        buildMaxTree(newSortedStringsLengths.begin() + i, newSortedStringsLengths.begin() + end, sortedStringsLengthsMaxTree[currentTable].begin());
      }

	    remainingAtLastDbUpdate = numberOfInputSequences - numberOfMarked;
	  } // else
	} // if ((numberOfInputSequences - numberOfMarked < remainingAtLastDbUpdate * 0.5) && (remainingAtLastDbUpdate > numberOfInputSequences * 0.05))
	querySequencePointer = lexicographicSequencePointers[lexicographicIndex];
	querySequenceLength = sequenceLengths[lexicographicIndex];
	queryHeaderPointer = headerPointers[lexicographicIndex];
	queryHeaderLength = headerLengths[lexicographicIndex];

	makeCluster();
	
	} // if (!markedGrey[lexicographicIndex])
      } // if (!marked[lexicographicIndex])


    } // for (Index i = 0; i < numberOfInputSequences; ++i)

    return EXIT_SUCCESS;
  }
} // main
