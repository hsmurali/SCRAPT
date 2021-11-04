#ifndef FASTA_HPP
#define FASTA_HPP

#include <string>
#include <istream>
#include <ostream>
#include <vector>

namespace sequence
{
  using std::string;
  using std::vector;
  using std::istream;
  using std::ostream;

  const int CHARS_PER_LINE = 60;
  
  struct FastaRecord
  {
    FastaRecord();
    FastaRecord(const string &);
    string header;
    string sequence;
  };

  istream &operator>>(istream &, FastaRecord &);
  ostream &operator<<(ostream &, const FastaRecord &);
  
  typedef vector<FastaRecord> Fasta;
  
  istream &operator>>(istream &, Fasta &);
  ostream &operator<<(ostream &, const Fasta &);
}

#endif
