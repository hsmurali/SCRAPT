#!/bin/bash
#SBATCH -J Benchmark_DNACLUST_Spatial_Soil # Job name
#SBATCH -o Benchmark_DNACLUST_Spatial_Soil.o%j # Name of output file
#SBATCH -e Benchmark_DNACLUST_Spatial_Soil.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=2-00:00:00
#SBATCH --qos=large
#SBATCH --mem=36gb
#SBATCH --cpus-per-task=16

data_path=/fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Earth_Microbiome/deduplicated.soil.seqs.fna
out_path=/fs/cbcb-lab/mpop/projects/SCRAPT/Experiments/Spatil_Soil/DNACLUST_Benchmarks/
prog_path=/fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/SCRAPT/dnaclust/dnaclust_linux_release3/dnaclust 

mkdir ${out_path}

sim=${1}

/usr/bin/time -v ${prog_path} ${data_path} -s ${sim} -t 8 --no-k-mer-filter  > ${out_path}dnaclust_${sim}.txt