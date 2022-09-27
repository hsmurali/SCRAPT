DNACLUST is a tool for clustering millions of short DNA sequences.

```DNACLUST``` uses ```gxx >= 9.3.0``` and run make after creating the conda environment.

To run ```DNACLUST```
```
bin/dnaclust -h 
Usage: DNACLUST [-h] --input-file VAR [--similarity VAR] [--predetermined-cluster-centers VAR] [--recruit-only] [--header] [--left-gaps-allowed] 
                [--k-mer-length VAR] [--approximate-filter] [--k-mer-filter] [--no-overlap] [--threads VAR] [--use-full-query-header] 
                [--mismatches VAR] [--assign-ambiguous] [--random-seed VAR] [--print-inverted-index]

The output is written to STDOUT.
Each line will contain the ids of the sequences in each cluster,and the first id of each line is the cluster representative.
Example: To cluster a set of 16S rRNA fragments at 0.98 similarity use:
bin/dnaclust -i file.fasta -l -s 0.98 > clusters 
You can optionally specify a k-mer length for the filter.The longer k-mers use more memory. Also the filter will be more specific with longer k-mers.The.default_value log_4(median length) should be good for most cases.


Optional arguments:
  -h, --help                         	shows help message and exits 
  -v, --version                      	prints version information and exits 
  -i, --input-file                   	A fasta file of the input sequences [required]
  -s, --similarity                   	set similarity between cluster center and cluster sequences [default: 0.99]
  -p, --predetermined-cluster-centers	file containing predetermined cluster centers [default: ""]
  -r, --recruit-only                 	when used with predetermined-cluster-centers option, only clusters the input sequences that are similar to the predetermined centers 
  -d, --header                       	output header line indicating run options 
  -l, --left-gaps-allowed            	allow for gaps on the left of shorter string in semi-global alignment 
  -k, --k-mer-length                 	length of k-mer for filtering [default: 0]
  --approximate-filter               	use faster approximate k-mer filter 
  --k-mer-filter                     	use k-mer filter 
  --no-overlap                       	cluster some of sequences such that the cluster centers are at distance at least two times the radius of the clusters 
  -t, --threads                      	Number of Threads [default: 1]
  -u, --use-full-query-header        	use the full query header instead of the first word 
  -m, --mismatches                   	number of mismatches allowed from cluster center [default: -1]
  -a, --assign-ambiguous             	assign ambiguous reads to clusters based on abundances of non-ambiguous reads 
  -e, --random-seed                  	Seed for random number generator [default: 0]
  --print-inverted-index             	Print mapping from sequence to each center
```

DNACLUST Copyright (C) 2010 Mohammadreza Ghodsi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
This program was developed by Mohammadreza Ghodsi and currently maintained by, Harihara Subrahmaniam Muralidharan. 
