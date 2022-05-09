#!/bin/bash
#SBATCH -J pre_processing # Job name
#SBATCH -o pre_processing.o.Arctic%j # Name of output file
#SBATCH -e pre_processing.e.Arctic%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --qos=workstation
#SBATCH --ntasks=12
#SBATCH --mem=48gb
#SBATCH --array=2-708


module load qiime/1.9.1
module load ea_utils
module load fastx_toolkit/0.0.14
module load fastqc

data_path=/fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Earth_Microbiome/Spatial_Soil/
out_dir=/fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Earth_Microbiome/split_libs_spatial_soil/
mkdir ${out_dir}

ls ${data_path} | grep "^ERR" > samples.txt
sample=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -n 1`
s=${sample}
echo ${s}
echo ${sample}

cd ${data_path}

split_libraries_fastq.py -i ${sample}/${sample}.fastq.gz -o ${out_dir}${s}/ -q 19 --phred_offset 33 --barcode_type 'not-barcoded' --sample_ids ${s}
