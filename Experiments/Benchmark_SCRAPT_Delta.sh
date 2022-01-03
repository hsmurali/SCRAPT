#!/bin/bash
#SBATCH -J Performance_Delta # Job name
#SBATCH -o Performance_Delta.o%j # Name of output file
#SBATCH -e Performance_Delta.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=large
#SBATCH --ntasks=8
#SBATCH --mem=36gb

module load anaconda
source activate /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/venvs/16S-Clustering/
alpha=${1}
delta=${2} #0.008
sim=0.95

	
data_path=/fs/cbcb-scratch/hsmurali/Iterative_Clustering/Datasets/Lupus-Microbiome-Unpublished/deduplicated.seqs.fna
counts_path=/fs/cbcb-scratch/hsmurali/Iterative_Clustering/Datasets/Lupus-Microbiome-Unpublished/Counts.dict
prog_path=/fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/SCRAPT/src/
outdir=/fs/cbcb-scratch/hsmurali/Iterative_Clustering/Experiments/Lupus-Microbiome-Unpublished/

mkdir ${outdir}
mkdir ${outdir}Parameterize_on_Delta/sim_${sim}/
mkdir ${outdir}Parameterize_on_Delta/sim_${sim}/alpha_${alpha}/

/usr/bin/time -v python ${prog_path}SCRAPT.py -f ${data_path} -c ${counts_path} -o ${outdir}Parameterize_on_Delta/sim_${sim}/alpha_${alpha}/delta_${delta}/ -s ${alpha} -a True -r ${sim} -m True -d ${delta}
