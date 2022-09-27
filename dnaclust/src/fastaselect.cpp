#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include "fasta.hpp"
#include "argparse.hpp"

template <class InputIterator>
std::string firstWord(InputIterator first, InputIterator last)
{
    return std::string(first, std::find_if(first, last, isspace));
}

int main(int argc, char *argv[])
{

    argparse::ArgumentParser program("fastaselect");
    program.add_description("Options\n");

    program.add_argument("-f", "--fasta-file")
        .help("A fasta file of the input sequences")
        .required();

    program.add_argument("-c", "--cluster-centers")
        .help("Write the sequences of all cluster centers to standard output in FASTA format.")
        .default_value(true).implicit_value(false);

    program.add_argument("-a", "--all-clusters")
        .help("Write the sequences of each cluster to a seperate FASTA file." 
            "The name of the files will be the given path and prefix,"
            "followed by cluster number andy a '.fasta' suffix.")
        .default_value(false).implicit_value(true);

    program.add_argument("-p", "--file-name-prefix")
        .help("Specify the path and prefix for the cluster FASTA file names.")
        .default_value(std::string{"/cluster-"});

    program.add_argument("--everything-except")
        .help("Write all sequences except cluster centers to the standard output in FASTA format.")
        .default_value(false)
        .implicit_value(true);

    program.add_argument("--cluster-sequences")
        .help("Write the cluster sequences (but not the cluster center sequence) to the standard output in FASTA format.")
        .default_value(false)
        .implicit_value(true);
    
    program.parse_args(argc, argv);

    bool cluster_centers = program.get<bool>("-c");
    bool all_clusters = program.get<bool>("-a");
    bool everything_except = program.get<bool>("--everything-except");
    bool cluster_sequences = program.get<bool>("--cluster-sequences");
    std::string fasta_file_name = program.get("-f");
    std::string clusters_file_name_prefix = program.get("-p");
    const std::string FASTA_SUFFIX = "fasta";

    std::ifstream fasta_file(fasta_file_name.c_str());
    sequence::Fasta sequences;
    fasta_file >> sequences;
    typedef std::map<std::string, size_t> StringIndexMap;
    StringIndexMap index_of_id;

    for (size_t i = 0; i < sequences.size(); ++i) 
    {
        const std::string &header = sequences[i].header;
        index_of_id[firstWord(header.begin(), header.end())] = i;
    }

    if (cluster_centers) 
    {
        std::string line;
        while (getline(std::cin, line)) 
        {
            std::string center_id = firstWord(line.begin(), line.end());
            std::cout << sequences[index_of_id[center_id]];
        }
    } 
    else if (everything_except) 
    {
        std::set<std::string> center_ids;
        {
            std::string line;
            while (getline(std::cin, line)) 
                center_ids.insert(firstWord(line.begin(), line.end()));
        }
        for (size_t i = 0; i < sequences.size(); ++i) 
        {
            const std::string &header = sequences[i].header;
            std::string id = firstWord(header.begin(), header.end());
            if (center_ids.find(id) == center_ids.end())
                std::cout << sequences[i];
        }
    } 
    else if (cluster_sequences) 
    {
        std::string line;
        while (getline(std::cin, line)) 
        {
            std::istringstream line_stream(line);
            std::string id;
            line_stream >> id;
            while (line_stream >> id)
                std::cout << sequences[index_of_id[id]];
        }
    }

    if (all_clusters) 
    {
        std::string line;
        int cluster_number = 0;
        while (getline(std::cin, line)) 
        {
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
