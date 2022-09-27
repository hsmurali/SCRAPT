#include "fasta.hpp"
#include <algorithm>
#include <ostream>
#include <stdexcept>

using std::min;
using std::endl;
using std::runtime_error;
const std::string WHITESPACE = " \n\r\t\f\v";
 
std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}
 
std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}

namespace sequence{

  FastaRecord::FastaRecord(const string &h): header(h)
  {
  }
  
  FastaRecord::FastaRecord()
  {
  }

  istream &operator>>(istream &inputStream, FastaRecord &fr)
  {
      fr.header.clear();
      fr.sequence.clear();

      if('>' != inputStream.get())
          throw runtime_error("Error in FASTA file format.");
      getline(inputStream, fr.header);
      while(inputStream.good() && '>' != inputStream.peek())
      {
          string line;
          getline(inputStream, line);
          fr.sequence += line;
      }
      return inputStream;
  }

  ostream &operator<<(ostream &outputStream, const FastaRecord &fr)
  {
      outputStream << '>' << fr.header;
      for(string::size_type i = 0; i < fr.sequence.length(); ++i)
      {
          if(0 == i % CHARS_PER_LINE)
              outputStream << '\n';
          outputStream << fr.sequence[i];
      }
      outputStream << endl;
      return outputStream;
  }
  
  istream &operator>>(istream &inputStream, Fasta &fa)
  {
      fa.clear();
      string line;
      while(getline(inputStream, line))
          if(not line.empty() and '>' == line[0])
              fa.push_back(FastaRecord(line.substr(1)));
          else
          {
              line = trim(line);
              if(not line.empty())
                  fa.rbegin()->sequence += line;
          }
      return inputStream;
  }

  ostream &operator<<(ostream &outputStream, const Fasta &fa)
  {
      for(Fasta::const_iterator i = fa.begin(); i != fa.end(); ++i)
          outputStream << *i;
      return outputStream;
  }  
}