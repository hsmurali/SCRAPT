#!/bin/bash
#SBATCH -J Mode_Against_Length # Job name
#SBATCH -o Mode_Against_Length.o # Name of output file
#SBATCH -e Mode_Against_Length.e # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=workstation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=48gb

module load anaconda
source activate /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/venvs/16S-Clustering/


data_path=/fs/cbcb-lab/mpop/mpop/Projects/QJRX/00.RawData/E20200314/E20200314_raw.fasta
counts_path=/fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Lupus-Microbiome-Published/Lengths.dict
prog_path=/fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/SCRAPT/src/
outdir=/fs/cbcb-scratch/hsmurali/Iterative_Clustering/E20200314/

rm -rf ${outdir}

/usr/bin/time -v python ${prog_path}SCRAPT.py -f ${data_path} -o ${outdir} -s 0.1 -a True -r 0.98 -m True -d 0.008 -t 8 -n 5
