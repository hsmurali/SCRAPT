#!/bin/bash
#SBATCH -J Iterative_Clustering_Delta # Job name
#SBATCH -o Iterative_Clustering_Delta.o%j # Name of output file
#SBATCH -e Iterative_Clustering_Delta.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --ntasks=8
#SBATCH --mem=36gb

module load anaconda
source activate /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/venvs/16S-Clustering/
alpha=${1}
delta=${2}

data_path=/fs/cbcb-scratch/hsmurali/16S-Clustering/data/Luo-Microbiome/Clusters_Implementation_2/deduplicated.seqs.fna
prog_path=/fs/cbcb-scratch/hsmurali/Iterative_Clustering/
outdir=/fs/cbcb-scratch/hsmurali/Iterative_Clustering/Experiments/

mkdir ${outdir}
mkdir ${outdir}/Vary_Delta/
mkdir ${outdir}/Vary_Delta/alpha_${alpha}/

adaptive=True

out=${outdir}/Vary_Delta/alpha_${alpha}/delta_${delta}

cd ${prog_path}
/usr/bin/time -v python SCRAPT.py -f ${data_path} -o ${out} -s ${alpha} -a True -d ${delta}

####Experiments to perform without modeshifting we performed by commenting out mode-shifting code. 