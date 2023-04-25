#PROFILE_FLAGS = -pg -fprofile-arcs -ftest-coverage -g -O1
CXXFLAGS = -Wall -std=c++17 -pedantic -O2 -static 
CXX = g++
#I/fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/venvs/SCRAPT/include

all: dnaclust fastaselect fastasort clean

dnaclust: dnaclust.o fasta.o
	$(CXX) $(CXXFLAGS) $(PROFILE_FLAGS) dnaclust.o fasta.o -o $@ -lpthread 

dnaclust.o: src/dnaclust.cpp src/search_include.cpp src/utility.hpp src/ternary_sort.hpp src/multi_dim.hpp
	$(CXX) -c $(CXXFLAGS) $(PROFILE_FLAGS) src/dnaclust.cpp

fastaselect: fastaselect.o fasta.o
	$(CXX) $(CXXFLAGS) $(PROFILE_FLAGS) fastaselect.o fasta.o -o $@

fastaselect.o: src/fastaselect.cpp
	$(CXX) -c $(CXXFLAGS) $(PROFILE_FLAGS) src/fastaselect.cpp

fastasort: fastasort.o fasta.o
	$(CXX) $(CXXFLAGS) $(PROFILE_FLAGS) $+ -o $@

fastasort.o: src/fastasort.cpp
	$(CXX) -c $(CXXFLAGS) $(PROFILE_FLAGS) $+ -o $@

fasta.o: src/fasta.cpp src/fasta.hpp
	$(CXX) -c $(CXXFLAGS) src/fasta.cpp

clean:
	rm -rf bin
	mkdir bin
	mv dnaclust bin/
	mv fastaselect bin/
	mv fastasort bin/
	rm *.o