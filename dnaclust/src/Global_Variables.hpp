typedef char CharType;
typedef const CharType * ConstCharPointer;
typedef std::vector<ConstCharPointer> ConstCharPointerVector;
typedef int PositionType;
typedef std::vector<PositionType> PositionVector;
typedef int SequenceNumber;
typedef intptr_t Index;
typedef int32_t CostType;
typedef std::vector<Index> IndexVector;
typedef std::vector<char> BoolVector;

const CharType endOfStringChar = '\0';
const CharType invalidChar = '^';
const CharType gapChar = '-';
const char HEADER_INDICATOR = '%';
const char CLUSTER_INDICATOR = '#';
const CostType MAX_COST = std::numeric_limits<CostType>::max() / 2;
const size_t MAX_THREADS = 60;
const size_t MAX_LENGTH = 4500;
const Index NULL_INDEX = -1;

float DEFAULT_SIMILARITY = 0.99;
float DEFAULT_MISMATCHES = -1.0;

//Inputs- From the user
bool outputMultipleAlignment = false;
bool leftGapsAllowed = false;
bool noKmerFilter = true;
bool kmerFilter = false;
bool approximateFilter2 = false;
bool noOverlap = false;
bool recruit_only = false;
bool useFullQueryHeader = false;
bool printInvertedIndex = false;
bool assignAmbiguous = false;
float similarity;
float mismatches;
int k_mer_length;
int random_seed = 0;

struct thread_args
{
    PositionType l;
    PositionType r;
    PositionType pos;
    Index intervalIndex;
    uint16_t tableIndex;
    int globalOffset;
};

struct Context {
    int id;
};

struct SearchResult
{
    SequenceNumber number;
    CostType cost;
};
typedef std::vector<SearchResult> SearchResultVector;

inline const CostType cost(const CharType c1, const CharType c2)
{
    return (c1 != c2) ? 1 : 0;
}

template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"{ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, " "));
    os<<"}";
    return os;
}