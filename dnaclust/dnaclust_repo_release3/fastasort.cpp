#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include "fasta.hpp"
#include <boost/program_options.hpp>		
namespace po = boost::program_options;


struct longer_sequence
{
  bool operator()(const sequence::FastaRecord &r1, const sequence::FastaRecord &r2)
  {
    return r1.sequence.length() > r2.sequence.length();
  }
};

int main(int argc, char *argv[])
{

  bool produce_help = false;
  bool do_random_shuffle = false;
  // bool lexicographic = false;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", po::bool_switch(&produce_help), "Produce help message.")
    ("random-shuffle,r", po::bool_switch(&do_random_shuffle), "Randomly shuffle sequences.")
    // ("lexicographic,l", po::bool_switch(&lexicographic), "Sort sequences in lexicographical order.")
    ;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  notify(vm);

  if (produce_help) {
    std::cerr << "Usage: fastasort [OPTIONS]\n"
	      << " The FASTA sequences are read from the STDIN.\n"
	      << " The sorted sequences are written to the STDOUT in FASTA format.\n"
	      << desc << '\n';
    exit(EXIT_FAILURE);
  }


  sequence::Fasta sequences;
  std::cin >> sequences;

  if (do_random_shuffle) {
    std::random_shuffle(sequences.begin(), sequences.end());
    std::cout << sequences;
  } else {
    // sort by length. the first sequence is the longest. 
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
