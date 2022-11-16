#include <sstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <limits>
#include <stack>
#include <fstream>
#include <stdint.h>
#include <pthread.h>

#include "argparse.hpp"
#include "ternary_sort.hpp"
#include "fasta.hpp"
#include "multi_dim.hpp"
#include "utility.hpp"
#include "kmer_utils.hpp"

enum BackPointer {LEFT, DOWN, DOWN_LEFT, ROOT};

typedef multi_dim::Matrix<KmerCount> CountMatrix;
typedef multi_dim::Matrix<CostType> DpTable;
typedef multi_dim::Matrix<BackPointer> BackPointerTable;

uint16_t running_threads = 0;
size_t num_threads = 1;
int threads_completed = 0;
pthread_mutex_t running_mutex;
pthread_t searchWorkers[MAX_THREADS];
struct Context searchWorkersContext[MAX_THREADS];
struct thread_args searchWorkersArgs[MAX_THREADS];
pthread_mutex_t searchWorkersMutexes[MAX_THREADS];
pthread_cond_t searchWorkersConds[MAX_THREADS];
pthread_cond_t searchWorkersCompleted;
pthread_mutex_t allThreadsMutex;

CountMatrix spectrumMatrix;
IndexVector spectrumSplitGuide, lexicographicIndexOfSortedStringsIndexes, lexicographicIndexOfSpectrumIndexes, spectrumSearchResults;
KmerSpectrum querySpectrum, queryLowSpectrum;
PositionVector spectrumRemainingGuide;
SpectrumIndexVector spectrumIndexes;

SequenceNumber numberOfInputSequences;
ConstCharPointerVector sortedStrings, lexicographicSequencePointers, headerPointers;
PositionVector headerLengths, sequenceLengths;
ConstCharPointer querySequencePointer, queryHeaderPointer;
PositionType querySequenceLength, queryHeaderLength;

std::vector< std::vector<CharType> > word(MAX_THREADS, std::vector<CharType>(MAX_LENGTH + 1));
std::vector< std::vector<PositionType> > wordMinLength(MAX_THREADS, std::vector<PositionType>(MAX_LENGTH + 1));
std::vector< std::vector<CostType> > minCost(MAX_THREADS, std::vector<CostType>(MAX_LENGTH + 1));
std::vector< std::vector<PositionType> > left(MAX_THREADS, std::vector<PositionType>(MAX_LENGTH + 1));
std::vector< std::vector<PositionType> > right(MAX_THREADS, std::vector<PositionType>(MAX_LENGTH + 1));
std::map< std::string, std::vector<std::string> > final_clusters;
std::map< std::string, std::vector<std::string> > inverted_index;
std::map< std::string, int > center_nonambig_counts;

SearchResultVector searchResults;
std::vector< SearchResultVector > threadSearchResults(MAX_THREADS);
std::vector< PositionVector > sortedStringsLengthsMinTree(MAX_THREADS);
std::vector< PositionVector > sortedStringsLengthsMaxTree(MAX_THREADS);
int resultsCounter[MAX_THREADS];

CostType radious;
DpTable *tables[MAX_THREADS];
BackPointerTable backTable(MAX_LENGTH + 1, MAX_LENGTH + 1);

void initialize()
{
  // initialize() implementation depends on whether we allow gaps on both ends.
    for(uint16_t i = 0; i < num_threads; ++i)
    {
        if (mismatches >= 0)
            left[i][0] = - static_cast<PositionType>(mismatches);
        else
            left[i][0] = - static_cast<PositionType>(floor(MAX_LENGTH * (1 - similarity)));
        right[i][0] = std::numeric_limits<PositionType>::max() / 2;
        std::fill(word[i].begin(), word[i].end(), invalidChar);

        tables[i] = new DpTable(MAX_LENGTH + 1, MAX_LENGTH + 1);
        (*tables[i])(0, 0) = 0;
        backTable(0, 0) = ROOT;

        for (size_t j = 1; j <  MAX_LENGTH + 1; ++j) 
        {
            (*tables[i])(0, j) = leftGapsAllowed ? 0 : (*tables[i])(0, j - 1) + cost(gapChar, invalidChar);    
            backTable(0, j) = LEFT;
        }
    }
} 

/*
  The implementation of search function is seperated in another file.
  The reason for this is to reuse the code to make two search functions:
  The first search function does not keep track of back pointers and would not find the alignment but is supposedly faster.
  The second search function will find the alignment.
*/

#include "search_include.cpp"

template <class InputIterator>
std::string firstWord(InputIterator first, InputIterator last)
{
    return std::string(first, std::find_if(first, last, isspace));
}

BoolVector marked;
BoolVector markedGrey;
int numberOfMarked;
int remainingAtLastDbUpdate;

sequence::Fasta inputFasta;
sequence::Fasta predeterminedCentersFasta;

template <class RandomAccessIterator>
RandomAccessIterator median(RandomAccessIterator first, RandomAccessIterator last)
{
    RandomAccessIterator mid = first + ((last - first) / 2);
    std::nth_element(first, mid, last);
    return mid;
}

Index twoMeans(const Index l, const Index r)
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

void clusterSort(const Index l, const Index r, const Index index)
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
    Index mid = twoMeans(l, r);
    cluster_information::sizeOfLeftCluster[index] = mid - l;
    clusterSort(l, mid, index + 1);
    clusterSort(mid, r, index + 2 * (mid - l));
}

IndexVector clusterSearchResults;

namespace layered 
{

    int layeredMaxMoreThresholds[MAXIMUM_K_MER_LENGTH];
    int layeredMinMoreThresholds[MAXIMUM_K_MER_LENGTH];

    int **allQueryCounts = 0;

    int layeredClusterSearch(Index index, int k)
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
                return layeredClusterSearch(index, k + 1);
            }
            for (Index i = cluster_information::l[index]; i < cluster_information::r[index]; ++i)
                clusterSearchResults.push_back(i);
            return cluster_information::r[index] - cluster_information::l[index];
        }
        return (layeredClusterSearch(index + 1, k) + layeredClusterSearch(index + 2 * cluster_information::sizeOfLeftCluster[index], k));
    }
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


void *threadWorker(void *context) 
{
    int id = (static_cast<struct Context *>(context))->id;
    pthread_mutex_lock(&searchWorkersMutexes[id]);
    while (1) 
    {
        pthread_cond_wait(&searchWorkersConds[id], &searchWorkersMutexes[id]);
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
    radious = static_cast<CostType>(floor((1 - similarity) * querySequenceLength));
    if (mismatches >= 0)
        radious = static_cast<CostType>(mismatches);
    if (noOverlap)
        radious *= 2;
    if (! noKmerFilter) {
        querySpectrum = countKmers(querySequencePointer, querySequenceLength);
    } 
    if (! noKmerFilter) 
    {   
        // Set allQueryCounts
        std::copy(querySpectrum.begin(), querySpectrum.end(), layered::allQueryCounts[k_mer_length - 1]);
        for (int k = 0; k < k_mer_length - 1; ++k)
            std::fill_n(layered::allQueryCounts[k], cluster_information::numberOfKmers[k], 0);
        for (int k = k_mer_length - 1; k > 0; --k)
            for (int j = 0; j < cluster_information::numberOfKmers[k]; ++j) {
                int j_prime = j / static_cast<int>(NUMBER_OF_NUCLEOTIDES);
                layered::allQueryCounts[k - 1][j_prime] += layered::allQueryCounts[k][j];
            }

        for (int k = 0; k < k_mer_length - 1; ++k)
            for (int j = 0; j < cluster_information::numberOfKmers[k]; ++j)
                if (layered::allQueryCounts[k][j] > MAX_COUNT)
                    layered::allQueryCounts[k][j] = MAX_COUNT;
    } 

    clusterSearchResults.clear();
    if (! noKmerFilter) 
    {
        {
            for (int i = 0; i < k_mer_length; ++i)
            if (radious == 0) 
            {
                layered::layeredMaxMoreThresholds[i] = 0;
                layered::layeredMinMoreThresholds[i] = 0;
            } 
            else 
            {
                layered::layeredMaxMoreThresholds[i] = std::max(radious, static_cast<CostType>(querySequenceLength * 0.04)) * (i + 1);
                layered::layeredMinMoreThresholds[i] = radious * (i + 1);
            }
        } 
        layered::layeredClusterSearch(0, 0);
    } // if (! noKmerFilter)
    IndexVector unmarkeds;
    if (! noKmerFilter) 
    {
        // Remove already marked sequences.
        unmarkeds.reserve(clusterSearchResults.size());
        for(IndexVector::const_iterator it = clusterSearchResults.begin(); it != clusterSearchResults.end(); ++it){
            Index lexicographicIndex = spectrumIndexes[*it].second;
            if (!marked[lexicographicIndex])
                unmarkeds.push_back(lexicographicIndex);
        }
        const Index numberOfUnmarkeds = unmarkeds.size();
        std::sort(unmarkeds.begin(), unmarkeds.end());
        sortedStrings.resize(numberOfUnmarkeds);
        PositionVector sortedStringsLengths(sortedStrings.size());
        for (Index i = 0; i < numberOfUnmarkeds; ++i) 
        {
            sortedStrings[i] = lexicographicSequencePointers[unmarkeds[i]];
            sortedStringsLengths[i] = sequenceLengths[unmarkeds[i]];
        }
        searchResults.clear();
        for(uint16_t i = 0; i < num_threads; ++i)
        {
            word[i][1] = invalidChar;
            wordMinLength[i][1] = MAX_LENGTH + 1;
        }    
        int interval_size = static_cast<int>(ceil(sortedStrings.size() / static_cast<double>(num_threads)));
        if (sortedStrings.size() < num_threads) {
            num_threads = 1;
            interval_size = sortedStrings.size();
        }

        running_threads = 0;
        for (unsigned int i = 0, currentTable = 0; i < sortedStrings.size(); i += interval_size, ++currentTable) 
        {
            unsigned int end = i + interval_size;
            if (i + interval_size > sortedStrings.size()) 
                end = sortedStrings.size();
            else 
                end = i + interval_size;
            sortedStringsLengthsMinTree[currentTable].resize(2 * (end - i));
            sortedStringsLengthsMaxTree[currentTable].resize(2 * (end - i));
            buildMinTree(sortedStringsLengths.begin() + i, sortedStringsLengths.begin() + end, sortedStringsLengthsMinTree[currentTable].begin());
            buildMaxTree(sortedStringsLengths.begin() + i, sortedStringsLengths.begin() + end, sortedStringsLengthsMaxTree[currentTable].begin());
        }
        if (numberOfUnmarkeds) 
        {
            running_threads = 0;
            threads_completed = 0;
            for (int i = 0; i < numberOfUnmarkeds; i += interval_size) 
            {
                ++running_threads;
                pthread_mutex_lock(&searchWorkersMutexes[running_threads-1]);
                struct thread_args *args = &searchWorkersArgs[running_threads-1];
                args->l = 0;
                if (i + interval_size > numberOfUnmarkeds) 
                    args->r = (numberOfUnmarkeds - 1) - i;
                else 
                    args->r = interval_size - 1;
                args->pos = 0;
                args->intervalIndex = 0;
                args->tableIndex = running_threads-1;
                args->globalOffset = i;
                pthread_cond_signal(&searchWorkersConds[running_threads-1]);
                pthread_mutex_unlock(&searchWorkersMutexes[running_threads-1]);
            }
            // Wait for all threads to complete.
            pthread_mutex_lock(&allThreadsMutex);
            while (threads_completed < running_threads)
                pthread_cond_wait(&searchWorkersCompleted, &allThreadsMutex);
            pthread_mutex_unlock(&allThreadsMutex);
            for(int i = 0; i < running_threads; i++)  {
                searchResults.insert(searchResults.end(), threadSearchResults[i].begin(), threadSearchResults[i].end());
                threadSearchResults[i].clear();
            }
        }
    } // if (! noKmerFilter)
    else 
    {
        for(uint16_t i = 0; i < num_threads; ++i){
            word[i][1] = invalidChar;
            wordMinLength[i][1] = MAX_LENGTH + 1;
        }
        searchResults.clear();
        int interval_size = static_cast<int>(ceil(sortedStrings.size() / static_cast<double>(num_threads)));
        if (sortedStrings.size() < num_threads) 
        {
            num_threads = 1;
            interval_size = sortedStrings.size();
        }
        running_threads = 0;
        threads_completed = 0;
        for (unsigned int i = 0; i < sortedStrings.size(); i += interval_size) 
        {
            ++running_threads;
            pthread_mutex_lock(&searchWorkersMutexes[running_threads - 1]);
            struct thread_args *args = &searchWorkersArgs[running_threads - 1];
            args->l = 0;
            if (i + interval_size > sortedStrings.size())
                args->r = (sortedStrings.size() - 1) - i;
            else
                args->r = interval_size - 1;
            args->pos = 0;
            args->intervalIndex = 0;
            args->tableIndex = running_threads-1;
            args->globalOffset = i;
            pthread_cond_signal(&searchWorkersConds[running_threads-1]);
            pthread_mutex_unlock(&searchWorkersMutexes[running_threads-1]);
        }
        pthread_mutex_lock(&allThreadsMutex);
        while (threads_completed < running_threads)
            pthread_cond_wait(&searchWorkersCompleted, &allThreadsMutex);
        pthread_mutex_unlock(&allThreadsMutex);
        for(int i = 0; i < running_threads; i++)  {
            searchResults.insert(searchResults.end(), threadSearchResults[i].begin(), threadSearchResults[i].end());
            threadSearchResults[i].clear();
        }
    }
    { 
        typedef std::vector<PositionType> CountsVector;
        // This vector will contain the maximum number of gaps that happen after position [i] in query string in ANY of the alignments to search results.
        CountsVector gapCounts(querySequenceLength);
        // This string will contain the query sequence with appropriate number of gaps inserted at each position to produce a multiple alignment of the cluster.
        std::string clusterCenterSequenceWithGaps;
        if (useFullQueryHeader) 
        {
            if (!assignAmbiguous)
                std::cout << queryHeaderPointer << '\t';
            else 
                if (!printInvertedIndex)
                    center_nonambig_counts[queryHeaderPointer] = 1;
        } 
        else 
        {
            if (!assignAmbiguous)
                std::cout << firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength) << '\t';
            else
                if (!printInvertedIndex)
                    center_nonambig_counts[firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength)] = 1;
        }
        for (SearchResultVector::const_iterator resultIterator = searchResults.begin(); resultIterator != searchResults.end(); ++resultIterator) 
        {
            Index lexicographicIndex;
            if (! noKmerFilter)
                lexicographicIndex = unmarkeds[resultIterator->number];
            else
                lexicographicIndex = lexicographicIndexOfSortedStringsIndexes[resultIterator->number];
            if (!assignAmbiguous)
                markedGrey[lexicographicIndex] = true;
            if ((!noOverlap) || resultIterator->cost < radious / 2) 
            {
                if (!marked[lexicographicIndex]) 
                {
                    if (!assignAmbiguous) {
                        marked[lexicographicIndex] = true;
                        ++numberOfMarked;
                    }
                    ConstCharPointer headerPointer = headerPointers[lexicographicIndex];
                    PositionType headerLength = headerLengths[lexicographicIndex];
                    if (useFullQueryHeader) {
                        if (!assignAmbiguous)
                            std::cout << headerPointer << '\t';
                        else {
                            if (printInvertedIndex)
                                std::cout << headerPointer << '\t' << queryHeaderPointer << '\n';
                            else {
                                inverted_index[headerPointer].push_back(queryHeaderPointer);
                                if (inverted_index[headerPointer].size() == 1)
                                    center_nonambig_counts[queryHeaderPointer] += 1;
                                else if (inverted_index[headerPointer].size() == 2)
                                    center_nonambig_counts[inverted_index[headerPointer][0]] -= 1;
                            }
                        }
                    } 
                    else {
                        if (!assignAmbiguous)
                            std::cout << firstWord(headerPointer, headerPointer + headerLength) << '\t';
                        else {
                            if (printInvertedIndex)
                                std::cout << firstWord(headerPointer, headerPointer + headerLength) << '\t' 
                                          << firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength) << '\n';
                            else {
                                inverted_index[firstWord(headerPointer, headerPointer + headerLength)].push_back(
                                               firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength));
                                if (inverted_index[firstWord(headerPointer, headerPointer + headerLength)].size() == 1)
                                    center_nonambig_counts[firstWord(queryHeaderPointer, queryHeaderPointer + queryHeaderLength)] += 1;
                                else if (inverted_index[firstWord(headerPointer, headerPointer + headerLength)].size() == 2) 
                                    center_nonambig_counts[inverted_index[firstWord(headerPointer, headerPointer + headerLength)][0]] -= 1;
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
    { 
        argparse::ArgumentParser program("DNACLUST");
        program.add_description("   The output is written to STDOUT.\n"
                        "Each line will contain the ids of the sequences in each cluster," 
                        "and the first id of each line is the cluster representative.\n"
                        "Example: To cluster a set of 16S rRNA fragments at 0.98 similarity use:\n"
                        "./dnaclust file.fasta -l -s 0.98 > clusters \n"
                        "You can optionally specify a k-mer length for the filter."
                        "The longer k-mers use more memory. "
                        "Also the filter will be more specific with longer k-mers."
                        "The.default_value log_4(median length) should be good for most cases.\n");

        program.add_argument("-i", "--input-file")
            .help("A fasta file of the input sequences")
            .required();

        program.add_argument("-s", "--similarity")
            .help("set similarity between cluster center and cluster sequences")
            .default_value(DEFAULT_SIMILARITY).scan<'f',float>();
        
        program.add_argument("-p","--predetermined-cluster-centers")
            .help("file containing predetermined cluster centers")
            .default_value(std::string{""});//.scan<'i',std::string>();

        program.add_argument("-r", "--recruit-only")
            .help("when used with predetermined-cluster-centers option, only clusters the input sequences that are similar to the predetermined centers")
            .default_value(recruit_only).implicit_value(true);

        program.add_argument("-d","--header")
            .help("output header line indicating run options")
            .default_value(false).implicit_value(true);
           
        program.add_argument("-l","--left-gaps-allowed")
            .help("allow for gaps on the left of shorter string in semi-global alignment")
            .default_value(leftGapsAllowed).implicit_value(!leftGapsAllowed);

        program.add_argument("-k","--k-mer-length")
            .help("length of k-mer for filtering")
            .default_value(0);

        program.add_argument("--approximate-filter")    
            .help("use faster approximate k-mer filter")
            .default_value(false).implicit_value(true);

        program.add_argument("--k-mer-filter")
            .help("use k-mer filter")
            .default_value(false).implicit_value(true);

        program.add_argument("--no-overlap")
            .help("cluster some of sequences such that the cluster centers are at distance at least two times the radius of the clusters")
            .default_value(false).implicit_value(true);

        program.add_argument("-t","--threads")
            .help("Number of Threads")
            .default_value(num_threads).scan<'i',int>();

        program.add_argument("-u","--use-full-query-header")
            .help("use the full query header instead of the first word")
            .default_value(false).implicit_value(true);

        program.add_argument("-m","--mismatches") 
            .help("number of mismatches allowed from cluster center")
            .default_value(DEFAULT_MISMATCHES).scan<'f',float>();

        program.add_argument("-a","--assign-ambiguous") 
            .help("assign ambiguous reads to clusters based on abundances of non-ambiguous reads")
            .default_value(assignAmbiguous).implicit_value(!assignAmbiguous);

        program.add_argument("-e","--random-seed")
            .help("Seed for random number generator")
            .default_value(random_seed).scan<'i',int>();

        program.add_argument("--print-inverted-index")
            .default_value(printInvertedIndex).implicit_value(!printInvertedIndex)
            .help("Print mapping from sequence to each center");
            
        program.parse_args(argc, argv);

        similarity = program.get<float>("-s");
        std::string inputFileName = program.get("-i");
        std::string predeterminedCentersFileName = program.get("-p");
        bool header = program.get<bool>("-d");
        recruit_only = program.get<bool>("-r");
        leftGapsAllowed = program.get<bool>("-l");
        k_mer_length = program.get<int>("-k");
        approximateFilter2 = program.get<bool>("approximate-filter");
        kmerFilter = program.get<bool>("k-mer-filter");
        noKmerFilter = !kmerFilter;    
        noOverlap = program.get<bool>("no-overlap");
        num_threads = program.get<int>("-t");
        useFullQueryHeader = program.get<bool>("-u");
        mismatches = program.get<float>("-m");
        assignAmbiguous = program.get<bool>("-a");
        random_seed = program.get<int>("-e");
        printInvertedIndex = program.get<bool>("print-inverted-index");
    
        if(inputFileName == "")
        {
            std::cerr << "Input File is Missing\n";
            exit(EXIT_FAILURE);
        }

        if(recruit_only == true && predeterminedCentersFileName == "")
        {
            std::cerr << "Recruit mode is enabled, but centers file missing\n\n";
            exit(EXIT_FAILURE);
        }

        std::ifstream inputFile(inputFileName.c_str());
        inputFile >> inputFasta;
        numberOfInputSequences = inputFasta.size();
        // Initialize all threads.
        for (size_t i = 0; i < MAX_THREADS; ++i) 
        {
            pthread_mutex_init(&searchWorkersMutexes[i], NULL);
            pthread_cond_init(&searchWorkersConds[i], NULL);
            //struct thread_args *args = &args_array[running_threads-1];
            searchWorkersContext[i].id = i;
            pthread_create(&searchWorkers[i], NULL, threadWorker, static_cast<void *>(&searchWorkersContext[i]));
        }
        pthread_mutex_init(&allThreadsMutex, NULL);
        pthread_cond_init(&searchWorkersCompleted, NULL);
        pthread_mutex_init(&running_mutex, NULL);
        
        if (recruit_only) 
        {
            std::ifstream predeterminedCentersFile(predeterminedCentersFileName.c_str());
            predeterminedCentersFile >> predeterminedCentersFasta;
        }
    
        if (header)
        { 
            std::cout << "Similarity\t\t\t" << similarity << std::endl;
            std::cout << "Recruit only\t\t\t" << recruit_only << std::endl;
            std::cout << "Left gaps allowed\t\t\t" << leftGapsAllowed << std::endl;
            std::cout << "K-mer length\t\t\t" << k_mer_length << std::endl;
            std::cout << "Approximate Filter\t\t\t" << approximateFilter2 << std::endl;
            std::cout << "kmer filter\t\t\t" << kmerFilter << std::endl;
            std::cout << "No kmer filter\t\t\t" << noKmerFilter << std::endl;
            std::cout << "No overlap\t\t\t" << noOverlap << std::endl;
            std::cout << "Num threads\t\t\t" << num_threads << std::endl;
            std::cout << "Use full query header\t\t\t" << useFullQueryHeader << std::endl;
            std::cout << "Mismatches\t\t\t" << mismatches << std::endl;
            std::cout << "Assign ambiguous\t\t\t" << assignAmbiguous << std::endl;
            std::cout << "Print inverted index\t\t\t" << printInvertedIndex << std::endl;
        }

        if (k_mer_length == 0) 
        {
            // Try to guess appropriate k-mer length based on median sequence length.
            typedef std::vector<size_t> size_vector;
            size_vector sizes;
            for (sequence::Fasta::const_iterator i = inputFasta.begin(); i != inputFasta.end(); ++i)
                sizes.push_back(i->sequence.length());
            size_t median_size = *median(sizes.begin(), sizes.end());
            k_mer_length = static_cast<int>(floor(log(static_cast<float>(median_size)) / log(static_cast<float>(NUMBER_OF_NUCLEOTIDES))));
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
            for (Index i = 0; i < numberOfInputSequences; ++i)
            {
                lexicographicSequencePointers[i] = inputFasta[i].sequence.c_str();
                pointerFastaIndexes[i].first = inputFasta[i].sequence.c_str();
                pointerFastaIndexes[i].second = i;
            }
            ternary_sort::ternarySort(lexicographicSequencePointers.data(), numberOfInputSequences);
            std::sort(pointerFastaIndexes.begin(), pointerFastaIndexes.end());
            PointerIndexVector pointerLexicographicIndexes(numberOfInputSequences);
            for (Index i = 0; i < numberOfInputSequences; ++i) 
            {
                pointerLexicographicIndexes[i].first = lexicographicSequencePointers[i];
                pointerLexicographicIndexes[i].second = i;
            }

            std::sort(pointerLexicographicIndexes.begin(), pointerLexicographicIndexes.end());
            IndexVector fastaIndexOfLexicographicIndexes(numberOfInputSequences);
            IndexVector lexicographicIndexOfFastaIndexes(numberOfInputSequences);
            for (Index i = 0; i < numberOfInputSequences; ++i) 
            {
                assert(pointerFastaIndexes[i].first == pointerLexicographicIndexes[i].first);
                fastaIndexOfLexicographicIndexes[pointerLexicographicIndexes[i].second] = pointerFastaIndexes[i].second;
                lexicographicIndexOfFastaIndexes[pointerFastaIndexes[i].second] = pointerLexicographicIndexes[i].second;
            }

            headerPointers.resize(numberOfInputSequences);
            headerLengths.resize(numberOfInputSequences);
            sequenceLengths.resize(numberOfInputSequences);

            for (Index i = 0; i < numberOfInputSequences; ++i) 
            {
                Index fastaIndex = fastaIndexOfLexicographicIndexes[i];
                headerPointers[i] = inputFasta[fastaIndex].header.c_str();
                headerLengths[i] = inputFasta[fastaIndex].header.length();
                sequenceLengths[i] = inputFasta[fastaIndex].sequence.length();
            }
        } // Calculate global vectors indexed by 'lexicographical index': lexicographicalSequencePointers, sequenceLengths, headerPointers, headerLengths

        typedef std::pair<PositionType, Index> LengthIndexPair;
        typedef std::vector<LengthIndexPair> LengthIndexVector;
        LengthIndexVector lengthLexicographicIndexes(numberOfInputSequences);
    
        for (Index i = 0; i < numberOfInputSequences; ++i) 
        {
            lengthLexicographicIndexes[i].first = sequenceLengths[i];
            lengthLexicographicIndexes[i].second = i;
        }

        std::sort(lengthLexicographicIndexes.begin(), lengthLexicographicIndexes.end(), std::greater<LengthIndexPair>());
        if (! noKmerFilter) 
        { 
            spectrumIndexes.resize(numberOfInputSequences);
            for (Index i = 0; i < numberOfInputSequences; ++i) 
            {
                spectrumIndexes[i].first = countKmers(lexicographicSequencePointers[i], sequenceLengths[i]);
                spectrumIndexes[i].second = i;
            }
            { 
                {
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
                {
                   clusterSort(0, numberOfInputSequences, 0);
                   remainingAtLastDbUpdate = numberOfInputSequences;
                }
            }
        } // if (! noKmerFilter)

        if (noKmerFilter) 
        {
            sortedStrings = lexicographicSequencePointers;
            lexicographicIndexOfSortedStringsIndexes.resize(sortedStrings.size());
            for (Index i = 0; i < static_cast<Index>(sortedStrings.size()); ++i)
                lexicographicIndexOfSortedStringsIndexes[i] = i;
        
            // Go through each thread and resize all the threads
            //int num_threads_to_use = num_threads;
            int interval_size = static_cast<int>(ceil(sortedStrings.size() / static_cast<double>(num_threads)));//sortedStrings.size() / num_threads_to_use;

            if (sortedStrings.size() < num_threads) 
            {
                num_threads = 1; //sortedStrings.size(); //1; //sortedStrings.size();
                interval_size = sortedStrings.size();//1;//sortedStrings.size();
            }
            running_threads = 0;
          
            for (unsigned int i = 0, currentTable = 0; i < sortedStrings.size(); i += interval_size, ++currentTable) 
            {
                unsigned int end = i + interval_size;
                if (i + interval_size > sortedStrings.size()) {
                    end = sortedStrings.size();
                }
                else {
                    end = i + interval_size;
                }
                sortedStringsLengthsMinTree[currentTable].resize(2 * (end - i));
                sortedStringsLengthsMaxTree[currentTable].resize(2 * (end - i));        
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
        for (sequence::Fasta::const_iterator predeterminedIterator = predeterminedCentersFasta.begin(); predeterminedIterator != predeterminedCentersFasta.end();  
            ++predeterminedIterator) 
        {
            querySequencePointer = predeterminedIterator->sequence.c_str();
            querySequenceLength = predeterminedIterator->sequence.length();
            queryHeaderPointer = predeterminedIterator->header.c_str();
            queryHeaderLength = predeterminedIterator->header.length();
            makeCluster();
        } // for

        if (recruit_only) 
        {
            // If we are assigning ambiguous reads based on abundance, handle it here.
            if (assignAmbiguous) 
            {
                // initialize random seed: */
                if (random_seed == 0) 
                {
                    srand(time(NULL));
                } 
                else {
                    srand(random_seed);
                }
                typedef std::map< std::string, std::vector<std::string> >::const_iterator MapIterator;
                for (MapIterator iter = inverted_index.begin(); iter != inverted_index.end(); iter++)
                {
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
                        final_clusters[bins[index]].push_back(iter->first);
                    }
                }
                for (MapIterator iter = final_clusters.begin(); iter != final_clusters.end(); iter++) 
                {
                    std::cout << iter->first << "\t";
                    typedef std::vector<std::string>::const_iterator ListIterator;
                    for (ListIterator seq_iter = iter->second.begin(); seq_iter != iter->second.end(); seq_iter++) 
                    {
                        std::cout << *seq_iter << "\t";
                    }
                    std::cout << "\n";
                }
            }
            return EXIT_SUCCESS;
        }
        for (Index i = 0; i < numberOfInputSequences; ++i) 
        {
            Index lexicographicIndex = lengthLexicographicIndexes[i].second;
            if (!marked[lexicographicIndex]) 
            {
                if (!markedGrey[lexicographicIndex]) 
                {
                    markedGrey[lexicographicIndex] = true;
                    marked[lexicographicIndex] = true;
                    ++numberOfMarked;
                    if ((numberOfInputSequences - numberOfMarked < remainingAtLastDbUpdate * 0.5) && 
                       (remainingAtLastDbUpdate > numberOfInputSequences * 0.05) && 
                       (numberOfInputSequences - numberOfMarked > 1)) 
                    {
                        if (! noKmerFilter) 
                        {
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
                            clusterSort(0, remainingAtLastDbUpdate, 0);
                        } // if (! noKmerFilter)
                        else 
                        { // No k-mer filter
                            ConstCharPointerVector newSortedStrings;
                            PositionVector newSortedStringsLengths;
                            IndexVector newLexicographicIndexOfSortedStringsIndexes;
                            for (Index i = 0; i < static_cast<Index>(sortedStrings.size()); ++i) 
                            {
                                Index lexicographicIndex = lexicographicIndexOfSortedStringsIndexes[i];
                                if (! marked[lexicographicIndex])
                                {
                                    newSortedStrings.push_back(lexicographicSequencePointers[lexicographicIndex]);
                                    newLexicographicIndexOfSortedStringsIndexes.push_back(lexicographicIndex);
                                    newSortedStringsLengths.push_back(sequenceLengths[lexicographicIndex]);
                                }
                            }
                            sortedStrings.swap(newSortedStrings);
                            lexicographicIndexOfSortedStringsIndexes.swap(newLexicographicIndexOfSortedStringsIndexes);
                            //size_t num_threads_to_use = num_threads;
                            int interval_size = static_cast<int>(ceil(sortedStrings.size() / static_cast<double>(num_threads)));//sortedStrings.size() / num_threads_to_use;
                            if (sortedStrings.size() < num_threads) {
                                num_threads = 1;//sortedStrings.size(); //1;
                                interval_size = sortedStrings.size();//1; //sortedStrings.size();
                            }
                            for (unsigned int i = 0, currentTable = 0; i < sortedStrings.size() ; i += interval_size, ++currentTable) 
                            {
                                unsigned int end = i + interval_size;
                                if (i + interval_size > sortedStrings.size()) {
                                    end = sortedStrings.size();
                                }
                                else {
                                    end = i + interval_size;
                                }
                                sortedStringsLengthsMinTree[currentTable].resize(2 * (end - i));
                                sortedStringsLengthsMaxTree[currentTable].resize(2 * (end - i));
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