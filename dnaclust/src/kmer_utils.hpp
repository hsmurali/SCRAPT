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