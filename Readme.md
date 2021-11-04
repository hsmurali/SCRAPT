```SCRAPT``` is tool that performs an iterative clustering of DNA sequences. SCRAPT: **S**ample **C**luster **R**ecruit **A**da**P**t and i**T**erate iteratively clusters by only clustering subsample and recruits more sequences into the clusters identified. SCRAPT identifies larger clusters much early on, only in a fraction of time taken by other conventional clustering methods. 

```
usage: SCRAPT.py [-h] -f FILEPATH -o OUTPUT_DIRECTORY [-s SAMPLING_RATE]
                 [-a ADAPTIVE] [-d DELTA] [-r SIMILARITY] [-m MAX_ITERATIONS]
                 [-k MIN_CLUSTER] [-t NUM_THREADS]

SCRAPT: Sampling Clustering Recruiting AdaPt and iTerate.SCRAPT is a tool to
cluster 16S genen sequences, using an iterative approach

optional arguments:
  -h, --help            show this help message and exit
  -f FILEPATH, --filepath FILEPATH
                        Path to the file containing 16S sequences. At the
                        moment we support only fatsa file containing the 16S
                        reads.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Location to write the outputs to
  -s SAMPLING_RATE, --Sampling_Rate SAMPLING_RATE
                        Initial Sampling Rate (between 0 and 100)
  -a ADAPTIVE, --adaptive ADAPTIVE
                        Flag to run SCRAPT in adaptve mode
  -d DELTA, --delta DELTA
                        Adjustment constant for the adaptive sampling
  -r SIMILARITY, --similarity SIMILARITY
                        Similairity to run clustering with
  -m MAX_ITERATIONS, --max_iterations MAX_ITERATIONS
                        Maximum number of iterations to run the iterative
                        clustering.
  -k MIN_CLUSTER, --min_cluster MIN_CLUSTER
                        Size of the smallest cluster to detect
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads```
                        
 SCRAPT uses DNACLUST internally to cluster and recruit sequences to cluster. 

### References
Ghodsi, M., Liu, B. & Pop, M. DNACLUST: accurate and efficient clustering of phylogenetic marker genes. BMC Bioinformatics 12, 271 (2011). https://doi.org/10.1186/1471-2105-12-271
