#!/bin/bash
#SBATCH -J Benchmark_CDHIT_Lupus_Microbiome # Job name
#SBATCH -o Benchmark_CDHIT_Lupus_Microbiome.o%j # Name of output file
#SBATCH -e Benchmark_CDHIT_Lupus_Microbiome.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=5-00:00:00
#SBATCH --qos=large
#SBATCH --mem=36gb
#SBATCH --cpus-per-task=16
#SBATCH --nodelist=tern00

module load cdhit

######Lupus-Microbiome:
	##seqs:   /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Lupus-Microbiome-Published/deduplicated.seqs.fna
	##outdir: /fs/cbcb-lab/mpop/projects/SCRAPT/Experiments/Lupus-Microbiome-Published/CDHIT_EST_Benchmarks_g_1/

######Earth-Microbiome:
	##seqs:   /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Earth_Microbiome/deduplicated.soil.seqs.fna
	##outdir: /fs/cbcb-lab/mpop/projects/SCRAPT/Experiments/Spatil_Soil/CDHIT_EST_Benchmarks_g_1/

data_path=${1}
out_path=${2}
sim=${3}

mkdir ${out_path}
mkdir ${out_path}CDHIT_EST_${sim}_g_1/
echo ${sim}

/usr/bin/time -v cd-hit-est -i ${data_path} -o ${out_path}CDHIT_EST_${sim}_g_1/CDHIT_EST.${sim}.clusters -c ${sim} -T 8 -M 3000 -n 10 -g 1 > Lupus_CDHIT_EST_${sim}_g_1.logs