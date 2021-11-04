#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include "fasta.hpp"
#include <boost/program_options.hpp>		
namespace po = boost::program_options;

template <class InputIterator>
std::string firstWord(InputIterator first, InputIterator last)
{
  return std::string(first, std::find_if(first, last, isspace));
}

int main(int argc, char *argv[])
{

  bool produce_help;
  bool cluster_centers;
  bool all_clusters;
  std::string fasta_file_name;
  std::string clusters_file_name_prefix;
  const std::string FASTA_SUFFIX = "fasta";

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", po::bool_switch(&produce_help), "Produce help message.")
    ("fasta-file,f", po::value<std::string>(&fasta_file_name), "REQUIRED - File containing sequences in FASTA format.")
    ("cluster-centers,c", po::bool_switch(&cluster_centers), "Write the sequences of all cluster centers to standard output in FASTA format.")
    ("all-clusters,a", po::bool_switch(&all_clusters), "Write the sequences of each cluster to a seperate FASTA file. The name of the files will be the given path and prefix, folloed by cluster number andy a \'.fasta\' suffix.")
    ("file-name-prefix,p", po::value<std::string>(&clusters_file_name_prefix)->default_value("./cluster-"), "Specify the path and prefix for the cluster FASTA file names.")
    ;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  notify(vm);

  if (produce_help || fasta_file_name.empty()) {
    std::cerr << "Usage: fastaselect FastaFile\n"
	      << " The ids are read from standard input. Each line corresponds to a cluster with the first id being the cluster center.\n"
	      << " The output sequences are written to standard output in FASTA format.\n"
	      << desc << '\n';
    exit(EXIT_FAILURE);
  }


  std::ifstream fasta_file(fasta_file_name.c_str());
  sequence::Fasta sequences;
  fasta_file >> sequences;

  typedef std::map<std::string, size_t> StringIndexMap;

  StringIndexMap index_of_id;

  for (size_t i = 0; i < sequences.size(); ++i) {
    const std::string &header = sequences[i].header;
    index_of_id[firstWord(header.begin(), header.end())] = i;
  }


  if (cluster_centers) {
    std::string line;
    while (getline(std::cin, line)) {
      std::string center_id = firstWord(line.begin(), line.end());
      std::cout << sequences[index_of_id[center_id]];
    }

  }

  if (all_clusters) {
    std::string line;
    int cluster_number = 0;
    while (getline(std::cin, line)) {
      std::ostringstream cluster_file_name_stream;
      cluster_file_name_stream << clusters_file_name_prefix << ++cluster_number << '.' << FASTA_SUFFIX;
      std::ofstream cluster_file(cluster_file_name_stream.str().c_str());

      std::istringstream line_stream(line);
      std::string id;
      while (line_stream >> id)
	cluster_file << sequences[index_of_id[id]];
      
    }
  }

  return EXIT_SUCCESS;
}
