#include "Global_Variables.hpp"

typedef int16_t KmerNumber;
typedef uint16_t KmerCount;
typedef std::vector<KmerCount> KmerSpectrum; // k-mer spectrum is a vector of integers counting how many of each k-mer there are in a sequence.
typedef std::vector<KmerCount> KmerCountVector;
typedef std::pair<KmerSpectrum, Index> SpectrumIndexPair;
typedef std::vector<SpectrumIndexPair> SpectrumIndexVector;
typedef const KmerCount * const ConstKmerCountArray;
typedef KmerCount * const KmerCountArray;

const KmerCount MAX_COUNT = std::numeric_limits<KmerCount>::max();
const KmerCount MIN_COUNT = 0;
const int MAXIMUM_K_MER_LENGTH = 6;

enum Nucleotides{A, C, G, T, NUMBER_OF_NUCLEOTIDES};
KmerNumber number_of_k_mers;
const KmerCount MAX_K_MER_COUNT = std::numeric_limits<KmerCount>::max();

template <typename BaseType, typename ExponentType>
BaseType power(const BaseType base, const ExponentType exponent)
{
    BaseType result = 1;
    for (ExponentType i = 0; i < exponent; ++i)
        result *= base;
    return result;
}

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
    switch(tolower(nucleotide)) 
    {
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

// Count the number of each k-mer in a sequence and return a vector of counts.
KmerSpectrum countKmers(ConstCharPointer sequence, PositionType length)
{
    KmerSpectrum spectrum(number_of_k_mers);
    if (length >= k_mer_length) 
    {
        int i = 0;
        int kMerNumber = 0;
        for (; i < k_mer_length; ++i) 
        {
            kMerNumber *= NUMBER_OF_NUCLEOTIDES;
            kMerNumber += numberFromNucleotide(sequence[i]);
        }
        ++spectrum[kMerNumber];
        for (; i < length; ++i) 
        {
            kMerNumber *= NUMBER_OF_NUCLEOTIDES;
            kMerNumber %= number_of_k_mers;
            kMerNumber += numberFromNucleotide(sequence[i]);
            if (spectrum[kMerNumber] < MAX_K_MER_COUNT)
                ++spectrum[kMerNumber];
        }    
    }
    return spectrum;
}

template <class Iterator1, class Iterator2>
Iterator2 buildMaxTree(Iterator1 first, Iterator1 last, Iterator2 intervalIterator)
{
    if (last > first) 
    {
        if (last == first + 1) {
            *intervalIterator = *first;
        } 
        else 
        {
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
    if (last > first) 
    {
        if (last == first + 1) 
        {
            *intervalIterator = *first;
        } 
        else 
        {
            Iterator1 mid = first + (last - 1 - first) / 2 + 1;
            *intervalIterator = std::min(*buildMinTree(first,  mid, intervalIterator + 1),
            *buildMinTree(mid,   last, intervalIterator + 2 * (mid - first)));
        }
    }
    return intervalIterator;
}

template <class RandomAccessIterator>
RandomAccessIterator median(RandomAccessIterator first, RandomAccessIterator last)
{
    RandomAccessIterator mid = first + ((last - first) / 2);
    std::nth_element(first, mid, last);
    return mid;
}

Index twoMeans(const Index l, const Index r, SpectrumIndexVector& spectrumIndexes)
{
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
    int clusterChanges = 0;


    KmerCountVector counts1;
    KmerCountVector counts2;

    counts1.reserve(r - l);
    counts2.reserve(r - l);

    while((clustersChanged > 0) && ((clusterChanges <= 1) || (clustersChanged > (r - l) / 10))) 
    {
        ++clusterChanges;
        clustersChanged = 0;
    
        for (Index i = l; i < r; ++i) 
        {
            int d1 = 0;
            int d2 = 0;

            ConstKmerCountArray spectrum = spectrumIndexes[i].first.data();
            for (KmerNumber j = 0; j < number_of_k_mers; ++j) {
                d1 += abs(static_cast<int>(center1[j]) - static_cast<int>(spectrum[j]));
                d2 += abs(static_cast<int>(center2[j]) - static_cast<int>(spectrum[j]));
            }
            if (d1 != d2) 
            {
                int newClusterNumber = (d1 < d2) ? 1 : 0;
                if (newClusterNumber != clusterNumber[i - l]) {
                    clusterNumber[i - l] = newClusterNumber;
                    ++clustersChanged;
                }
            }
        } 
        if (clustersChanged) 
        { // Find new centers.
            for (KmerNumber j = 0; j < number_of_k_mers; ++j) 
            {
                counts1.clear();
                counts2.clear();
                for (Index i = l; i < r; ++i) {
                    if (clusterNumber[i - l] == 1) {
                        counts1.push_back(spectrumIndexes[i].first[j]);
                    } 
                    else {
                        counts2.push_back(spectrumIndexes[i].first[j]);
                    } 
                }    
                center1[j] = *median(counts1.begin(), counts1.end());
                center2[j] = *median(counts2.begin(), counts2.end());
            } 
        } 
    } 

    { 
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
                std::swap(clusterNumber[i - l], clusterNumber[j - l]);
                ++i;
                --j;
            } 
        } 
    } 
    Index result = l;
    while (clusterNumber[result - l] == 1)
        ++result;
    for (Index i = result; i < r; ++i)
        assert(clusterNumber[i - l] == 0);
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

void clusterSort(const Index l, const Index r, const Index index, SpectrumIndexVector& spectrumIndexes)
{
    {
        for (int k = 0; k < k_mer_length; ++k) {
            std::fill_n(&cluster_information::allMinCounts[k][index * cluster_information::numberOfKmers[k]], cluster_information::numberOfKmers[k], MAX_COUNT);
            std::fill_n(&cluster_information::allMaxCounts[k][index * cluster_information::numberOfKmers[k]], cluster_information::numberOfKmers[k], MIN_COUNT);
        }
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
        } 
        cluster_information::l[index] = l;
        cluster_information::r[index] = r;
    }
    bool allEqual = true;
    {
        for (KmerNumber j = 0; j < number_of_k_mers; ++j) 
        {
            KmerCount &minCount = cluster_information::allMinCounts[k_mer_length - 1][index * cluster_information::numberOfKmers[k_mer_length - 1] + j];
            KmerCount &maxCount = cluster_information::allMaxCounts[k_mer_length - 1][index * cluster_information::numberOfKmers[k_mer_length - 1] + j];
            if (minCount < maxCount)    
                allEqual = false;
        }
    }
    if (allEqual) {
        cluster_information::sizeOfLeftCluster[index] = NULL_INDEX;
        return;
    }
    Index mid = twoMeans(l, r, spectrumIndexes);
    cluster_information::sizeOfLeftCluster[index] = mid - l;
    clusterSort(l, mid, index + 1, spectrumIndexes);
    clusterSort(mid, r, index + 2 * (mid - l), spectrumIndexes);
}
