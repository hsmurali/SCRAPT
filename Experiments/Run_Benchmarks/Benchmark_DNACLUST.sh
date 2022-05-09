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
#SBATCH --nodelist=tern00

######Lupus-Microbiome:
	##seqs:   /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Lupus-Microbiome-Published/deduplicated.seqs.fna
	##outdir: /fs/cbcb-lab/mpop/projects/SCRAPT/Experiments/Lupus-Microbiome-Published/DNACLUST_Benchmarks/

######Earth-Microbiome:
	##seqs:   /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Earth_Microbiome/deduplicated.soil.seqs.fna
	##outdir: /fs/cbcb-lab/mpop/projects/SCRAPT/Experiments/Spatil_Soil/DNACLUST_Benchmarks/

prog_path=/fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/SCRAPT/dnaclust/dnaclust_linux_release3/dnaclust 
data_path=${1}
out_path=${2}
sim=${3}

mkdir ${out_path}

/usr/bin/time -v ${prog_path} ${data_path} -s ${sim} -t 8 --no-k-mer-filter  > ${out_path}dnaclust_${sim}.txt