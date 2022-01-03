#!/bin/bash
#SBATCH -J Benchmark_DNACLUST # Job name
#SBATCH -o Benchmark_DNACLUST.o%j # Name of output file
#SBATCH -e Benchmark_DNACLUST.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --ntasks=8

data_path=/fs/cbcb-scratch/hsmurali/Iterative_Clustering/Datasets/Lupus-Microbiome-Unpublished/deduplicated.seqs.fna
out_path=/fs/cbcb-scratch/hsmurali/Iterative_Clustering/Experiments/Lupus-Microbiome-Unpublished/DNACLUST_Benchmarks/
mkdir ${out_path}
s
sim=${1}

/usr/bin/time -v /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/SCRAPT/dnaclust/dnaclust_linux_release3/dnaclust  ${data_path} -s ${sim} -t 8 -k 0 --no-k-mer-filter  > ${out_path}dnaclust_${sim}.txt