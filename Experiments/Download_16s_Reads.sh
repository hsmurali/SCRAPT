#!/bin/bash
#SBATCH -J Download # Job name
#SBATCH -o Download.o%j # Name of output file
#SBATCH -e Download.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-1:00:00
#SBATCH --qos=large
#SBATCH --mem=36gb
#SBATCH --array=2-130

module load ncbi-sra
cd /fs/cbcb-scratch/hsmurali/Lupus-16S-Dataset_130

f=`head -n ${SLURM_ARRAY_TASK_ID} SRR_Acc_List_2.txt | tail -n 1`

mkdir ${f}
cd ${f}
fastq-dump --split-files ${f}
