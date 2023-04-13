```SCRAPT```: **S**ample **C**luster **R**ecruit **A**da**P**t and i**T**erate is a tool that performs an iterative clustering of DNA sequences.  We propose an iterative sampling-based 16S rRNA sequence clustering approach that targets the largest clusters in the data set, allowing users to stop the clustering process when sufficient clusters are available for the specific analysis being targeted. We describe a probabilistic analysis of the iterative clustering process that supports the intuition that the clustering process identifies the larger clusters in the data set first. Using a real data set of 16S rRNA gene sequence data, we show that the iterative algorithm, coupled with an adaptive sampling process and a mode-shifting strategy for identifying cluster representatives, substantially speeds up the clustering process while being effective at capturing the large clusters in the dataset. The experiments also shows ```SCRAPT``` is able to produce OTUs which are less fragmented than popular tools: ```UCLUST```, ```CD-HIT``` and  ```DNACLUST```.

## Software Requirements:

```SCRAPT``` is written in Python 3 and uses the following python packages. 
1. python >= 3.7.10
2. seqkit >= 0.16.0
3. pandas >= 1.2.4
4. numpy >= 1.20.2
5. gxx >= 9.3.0

An Environment.yml is also available and a conda environment can be created as,
```
conda env create -f Environment.yml --prefix <path-to-install>
```

On installing the conda environment, activate it as
```
conda activate <path-to-install>
```

```SCRAPT``` uses ```DNACLUST``` as the default clustering kernel.  To install ```DNACLUST```, run the following commands. (Please make sure the virtual environment is active before compiling ```DNACLUST```.
```
cd SCRAPT/dnaclust
make
```
On installing ```DNACLUST```, run ```SCRAPT``` as, 

```
usage: SCRAPT.py [-h] -f FILEPATH -o OUTPUT_DIRECTORY [-s SAMPLING_RATE]
                 [-a ADAPTIVE] [-d DELTA] [-r SIMILARITY] [-n MAX_ITERATIONS]
                 [-k MIN_CLUSTER] [-t NUM_THREADS] [-c COUNTS_DICT]
                 [-m MODE_SHIFT] [-p PROB] [-b REALIZATIONS]

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
  -p PROB, --prob PROB  Confidence threshold. [DEFAULT = 0.9]
  -b REALIZATIONS, --realizations REALIZATIONS
                        Number of realizations to estimate confidence bounds.
                        [DEFAULT = 1000]
```

[![DOI](https://zenodo.org/badge/424442689.svg)](https://zenodo.org/badge/latestdoi/424442689)

## References
Ghodsi, M., Liu, B. & Pop, M. DNACLUST: accurate and efficient clustering of phylogenetic marker genes. BMC Bioinformatics 12, 271 (2011). https://doi.org/10.1186/1471-2105-12-271

If you use ```SCRAPT``` for your research, please cite the article published in Nucleic Acids Research.
```
@article{10.1093/nar/gkad158,
    author = {Luan, Tu and Muralidharan, Harihara Subrahmaniam and Alshehri, Marwan and Mittra, Ipsa and Pop, Mihai},
    title = "{SCRAPT: an iterative algorithm for clustering large 16S rRNA gene data sets}",
    journal = {Nucleic Acids Research},
    year = {2023},
    month = {03},
    abstract = "{16S rRNA gene sequence clustering is an important tool in characterizing the diversity of microbial communities. As 16S rRNA gene data sets are growing in size, existing sequence clustering algorithms increasingly become an analytical bottleneck. Part of this bottleneck is due to the substantial computational cost expended on small clusters and singleton sequences. We propose an iterative sampling-based 16S rRNA gene sequence clustering approach that targets the largest clusters in the data set, allowing users to stop the clustering process when sufficient clusters are available for the specific analysis being targeted. We describe a probabilistic analysis of the iterative clustering process that supports the intuition that the clustering process identifies the larger clusters in the data set first. Using real data sets of 16S rRNA gene sequences, we show that the iterative algorithm, coupled with an adaptive sampling process and a mode-shifting strategy for identifying cluster representatives, substantially speeds up the clustering process while being effective at capturing the large clusters in the data set. The experiments also show that SCRAPT (Sample, Cluster, Recruit, AdaPt and iTerate) is able to produce operational taxonomic units that are less fragmented than popular tools: UCLUST, CD-HIT and DNACLUST. The algorithm is implemented in the open-source package SCRAPT. The source code used to generate the results presented in this paper is available at https://github.com/hsmurali/SCRAPT.}",
    issn = {0305-1048},
    doi = {10.1093/nar/gkad158},
    url = {https://doi.org/10.1093/nar/gkad158},
    note = {gkad158},
    eprint = {https://academic.oup.com/nar/advance-article-pdf/doi/10.1093/nar/gkad158/49515274/gkad158.pdf},
}
```
