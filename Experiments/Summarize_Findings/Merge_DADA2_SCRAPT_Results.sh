#!/bin/bash
#SBATCH -J Merge_DADA2_SCRAPT # Job name
#SBATCH -o Merge_DADA2_SCRAPT.o%j # Name of output file
#SBATCH -e Merge_DADA2_SCRAPT.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb

module load anaconda
source activate /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/venvs/SCRAPT

scrapt_sim=${1}
alpha=${2}
dada2_sim=${3}

python /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/SCRAPT/Experiments/Summarize_Findings/Merge_DADA2_SCRAPT_Results.py ${scrapt_sim} ${alpha} ${dada2_sim}