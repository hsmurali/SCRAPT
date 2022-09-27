#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include "fasta.hpp"
#include "argparse.hpp"		

struct longer_sequence
{
    bool operator()(const sequence::FastaRecord &r1, const sequence::FastaRecord &r2)
    {
        return r1.sequence.length() > r2.sequence.length();
    }
};

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program("fastasort");
    program.add_description("Options\n");
    program.add_argument("-r","--random-shuffle")
        .help("Randomly shuffle sequences.")
        .default_value(false)
        .implicit_value(true);
    program.parse_args(argc, argv);
    bool do_random_shuffle = program.get<bool>("-r");

    sequence::Fasta sequences;
    std::cin >> sequences;

    if (do_random_shuffle) 
    {
        std::random_shuffle(sequences.begin(), sequences.end());
        std::cout << sequences;
    } 
    else 
    {
        typedef std::pair<size_t, size_t> length_index_pair;
        std::vector<length_index_pair> length_indexes;
        for (size_t i = 0; i < sequences.size(); ++i)
            length_indexes.push_back(std::make_pair(sequences[i].sequence.length(), i));
        std::sort(length_indexes.begin(), length_indexes.end(), std::greater<length_index_pair>());
        for (size_t i = 0; i < length_indexes.size(); ++i)
            std::cout << sequences[length_indexes[i].second];
    }

    return EXIT_SUCCESS;
}
