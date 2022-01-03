#!/bin/bash
#SBATCH -J pre_processing # Job name
#SBATCH -o pre_processing.o # Name of output file
#SBATCH -e pre_processing.e # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --qos=workstation
#SBATCH --mem=48gb

module load qiime/1.9.1
module load ea_utils
module load fastx_toolkit/0.0.14
module load fastqc

cd /fs/cbcb-scratch/hsmurali/Lupus-16S-Dataset

# JOIN PAIRED END READS
multiple_join_paired_ends.py -i samples -o joined_fastq --read1_indicator '_1' --read2_indicator '_2'

for file in joined_fastq/*_1; do
	mv $file ${file/_1/}
done

echo -e 'SampleID\tMergedReadCount' > joined_fastq_read_counts.txt
for file in $(ls joined_fastq); do
	count=`cat ./joined_fastq/${file}/fastqjoin.join.fastq | wc -l`
	echo -e ${file}'\t'$(( count / 4)) >> joined_fastq_read_counts.txt
done

# Remove unjoined reads so they don't mess up split librariess
rm joined_fastq/*/*un*

echo "split_libraries_fastq:phred_quality_threshold 19" > params.txt
echo "split_libraries_fastq:store_qual_scores True" >> params.txt
echo "split_libraries_fastq:phred_offset 33" >> params.txt
echo "split_libraries_fastq:barcode_type not-barcoded" >> params.txt

 # Merge into one clean fasta files
multiple_split_libraries_fastq.py -i joined_fastq -o split_libs -p params.txt --read_indicator fastqjoin.join --demultiplexing_method sampleid_by_file --include_input_dir_path --remove_filepath_in_name

# convert qual file to two lines per seq
fasta_formatter -i split_libs/seqs.qual -o split_libs/tmp.qual -w 0

# Remove extra brackets
cat ./split_libs/tmp.qual | sed 's/^\[//g' | sed 's/\]$//g' > ./split_libs/seqs_formatted.qual
rm ./split_libs/tmp.qual

convert_fastaqual_fastq.py -f ./split_libs/seqs.fna -q ./split_libs/seqs_formatted.qual -o ./split_libs/

fastqc ./split_libs/seqs.fastq -o ./split_libs
