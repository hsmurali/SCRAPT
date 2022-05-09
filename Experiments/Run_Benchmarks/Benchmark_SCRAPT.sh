#!/bin/bash
#SBATCH -J Benchmark_SCRAPT_Lupus_Unpub # Job name
#SBATCH -o Benchmark_SCRAPT_Lupus_Unpub.o%j # Name of output file
#SBATCH -e Benchmark_SCRAPT_Lupus_Unpub.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=12:00:00
#SBATCH --qos=large 
#SBATCH --mem=36gb
#SBATCH --cpus-per-task=16
#SBATCH --nodelist=tern00

module load anaconda
source activate /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/venvs/16S-Clustering/

######Lupus-Microbiome:
	##seqs:   /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Lupus-Microbiome-Published/deduplicated.seqs.fna
	##counts: /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Lupus-Microbiome-Published/Counts.dict
	##outdir: /fs/cbcb-lab/mpop/projects/SCRAPT/Experiments/Lupus_Microbiome_MT/

######Earth-Microbiome:
	##seqs:   /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Earth_Microbiome/deduplicated.soil.seqs.fna
	##counts: /fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Earth_Microbiome/Counts.soil.dict
	##outdir: /fs/cbcb-lab/mpop/projects/SCRAPT/Experiments/Spatil_Soil/

data_path=${1}
counts_path=${2}
outdir=${3}
alpha=${4}
sim=${5}
delta=0.008

prog_path=/fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/SCRAPT/src/

mkdir ${outdir}
mkdir ${outdir}Adaptive_With_Modeshifting/
mkdir ${outdir}Adaptive_Without_Modeshifting/
mkdir ${outdir}Fixed_With_Modeshifting/
mkdir ${outdir}Fixed_Without_Modeshifting/

mkdir ${outdir}Adaptive_With_Modeshifting/sim_${sim}/
mkdir ${outdir}Adaptive_Without_Modeshifting/sim_${sim}/
mkdir ${outdir}Fixed_With_Modeshifting/sim_${sim}/
mkdir ${outdir}Fixed_Without_Modeshifting/sim_${sim}/

/usr/bin/time -v python ${prog_path}SCRAPT.py -f ${data_path} -c ${counts_path} -o ${outdir}Adaptive_With_Modeshifting/sim_${sim}/alpha_${alpha}/ -s ${alpha} -a True -r ${sim} -m True -t 8
/usr/bin/time -v python ${prog_path}SCRAPT.py -f ${data_path} -c ${counts_path} -o ${outdir}Adaptive_Without_Modeshifting/sim_${sim}/alpha_${alpha}/ -s ${alpha} -a True -r ${sim} -m False -t 8
/usr/bin/time -v python ${prog_path}SCRAPT.py -f ${data_path} -c ${counts_path} -o ${outdir}Fixed_With_Modeshifting/sim_${sim}/alpha_${alpha}/ -s ${alpha} -a False -r ${sim} -m True -t 8
/usr/bin/time -v python ${prog_path}SCRAPT.py -f ${data_path} -c ${counts_path} -o ${outdir}Fixed_Without_Modeshifting/sim_${sim}/alpha_${alpha}/ -s ${alpha} -a False -r ${sim} -m False -t 8
