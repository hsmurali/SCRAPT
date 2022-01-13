```SCRAPT``` is tool that performs an iterative clustering of DNA sequences. SCRAPT: **S**ample **C**luster **R**ecruit **A**da**P**t and i**T**erate iteratively clusters by only clustering subsample and recruits more sequences into the clusters identified. SCRAPT identifies larger clusters much early on, only in a fraction of time taken by other conventional clustering methods. 

## Software Requirements:

```SCRAPT``` is written in Python 3 and uses the following python packages. 
1. python >= 3.7.10
2. seqkit >= 0.16.0
3. pandas >= 1.2.4
4. numpy >= 1.20.2
In addition to this, ```SCRAPT``` uses ```DNACLUST``` as an internal engine for clustering and searching. 

An Environment.yml is also available and a conda environment can be created as,
```
conda env create -f Environment.yml --prefix <path-to-install>
```
On installing the conda environment activate it through the following command,
```
conda activate <path-to-install>
```


```
usage: SCRAPT.py [-h] -f FILEPATH -o OUTPUT_DIRECTORY [-s SAMPLING_RATE]
                 [-a ADAPTIVE] [-d DELTA] [-r SIMILARITY] [-n MAX_ITERATIONS]
                 [-k MIN_CLUSTER] [-t NUM_THREADS] [-c COUNTS_DICT]
                 [-m MODE_SHIFT]

SCRAPT: Sampling Clustering Recruiting AdaPt and iTerate.SCRAPT is a tool to
cluster phylogenetic marker gene sequences, using an iterative approach.

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  -f FILEPATH, --filepath FILEPATH
                        Path to the file containing 16S sequences. At the
                        moment we support only fatsa file containing the 16S
                        reads.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Location to write the outputs to

optional named arguments:
  -s SAMPLING_RATE, --Sampling_Rate SAMPLING_RATE
                        Initial Sampling Rate (between 0 and 100). [DEFAULT =
                        0.1]
  -a ADAPTIVE, --adaptive ADAPTIVE
                        Flag to run SCRAPT in adaptve mode. [DEFAULT = True]
  -d DELTA, --delta DELTA
                        Adjustment constant for the adaptive sampling.
                        [DEFAULT = 0.008]
  -r SIMILARITY, --similarity SIMILARITY
                        Similairity threshold for clustering. [DEFAULT = 0.99]
  -n MAX_ITERATIONS, --max_iterations MAX_ITERATIONS
                        Maximum number of iterations to run the iterative
                        clustering. [DEFAULT = 50]
  -k MIN_CLUSTER, --min_cluster MIN_CLUSTER
                        Size of the smallest cluster to detect. [DEFAULT = 50]
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads. [DEFAULT = 8]
  -c COUNTS_DICT, --counts_dict COUNTS_DICT
                        Path to the dictionary of counts. If it is not
                        provided SCRAPT deduplicates the sequences.
  -m MODE_SHIFT, --mode_shift MODE_SHIFT
                        Perform Modeshifting. [DEFAULT = True]
```
                        
SCRAPT uses DNACLUST internally to cluster and recruit sequences to cluster. 

## References
Ghodsi, M., Liu, B. & Pop, M. DNACLUST: accurate and efficient clustering of phylogenetic marker genes. BMC Bioinformatics 12, 271 (2011). https://doi.org/10.1186/1471-2105-12-271
