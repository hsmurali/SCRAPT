#!/bin/bash
#SBATCH -J UCLUST-98-soil_8 # Job name
#SBATCH -o UCLUST-98-soil_8.o # Name of Output File
#SBATCH -e UCLUST-98-soil_8.e # Name of Error File
#SBATCH --mail-user=tluan@terpmail.umd.edu # Email for Job Info
#SBATCH --mail-type=all
#SBATCH --time=6-00:00:00
#SBATCH --qos=large
#SBATCH --ntasks=8
#SBATCH --mem=128G
#SBATCH --nodelist=tern00

module load usearch

######Lupus-Microbiome:
	##seqs:   /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Lupus-Microbiome-Published/deduplicated.seqs.fna
	##outdir: /fs/cbcb-lab/mpop/projects/SCRAPT/Experiments/Lupus-Microbiome-Published/UCLUST_Benchmarks/

######Earth-Microbiome:
	##seqs:   /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Earth_Microbiome/deduplicated.soil.seqs.fna
	##outdir: /fs/cbcb-lab/mpop/projects/SCRAPT/Experiments/Spatil_Soil/UCLUST_Benchmarks/

data_path=${1}
out_path=${2}
sim=${3}

mkdir ${out_path}
mkdir ${out_path}cluster${sim}/

/usr/bin/time -v usearch9.2.64  -cluster_fast ${data_path}  -threads 8 -id ${sim} -fulldp -clusters {out_path}cluster${sim}/c_ >{out_path}uclust-out-${sim}