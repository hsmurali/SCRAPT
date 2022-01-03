#PROFILE_FLAGS = -pg -fprofile-arcs -ftest-coverage -g -O1
CXXFLAGS = -Wall -pedantic -std=c++98 -O2 -static
CXX = g++ 
#-DCHRIS_DEBUG

all: dnaclust fastaselect fastasort

dnaclust: dnaclust.o fasta.o
	$(CXX) $(CXXFLAGS) $(PROFILE_FLAGS) dnaclust.o fasta.o -lboost_program_options -o $@ -lpthread

dnaclust.o: dnaclust.cpp search_include.cpp utility.hpp ternary_sort.hpp multi_dim.hpp
	$(CXX) -c $(CXXFLAGS) $(PROFILE_FLAGS) dnaclust.cpp

fastaselect: fastaselect.o fasta.o
	$(CXX) $(CXXFLAGS) $(PROFILE_FLAGS) fastaselect.o fasta.o -lboost_program_options -o $@

fastaselect.o: fastaselect.cpp
	$(CXX) -c $(CXXFLAGS) $(PROFILE_FLAGS) fastaselect.cpp

fastasort: fastasort.o fasta.o
	$(CXX) $(CXXFLAGS) $(PROFILE_FLAGS) $+ -lboost_program_options -o $@

fastasort.o: fastasort.cpp
	$(CXX) -c $(CXXFLAGS) $(PROFILE_FLAGS) $+ -o $@

fasta.o: fasta.cpp fasta.hpp
	$(CXX) -c $(CXXFLAGS) fasta.cpp

clean:
	rm -f dnaclust.o dnaclust fasta.o fastaselect.o fastaselect fastasort fastasort.o
