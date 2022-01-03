#!/bin/bash
similarity=0.98
threads=1

print_help()
{
    fold --spaces <<EOF
Usage: dnaclust-ref [OPTIONS...]
DNACLUST helper script to cluster sequences using a reference database.

  -c CENTERS         Fasta file of cluster centers/references.
  -d                 After clustering with reference database, perform de novo clustering.
  -r SIMILARITY      Set similarity between cluster center and cluster
                     sequences (default=0.98)
  -t THREADS         Set the number of threads to use
  -i INPUT_FILE      Fasta file of sequences to be clustered.
  -v                 Print verbose messages to standard error
  -h                 Give this help list

The sequences to be clustered are read from the STDIN. The cluster centers are written to STDOUT. Messages are written to STDERR.
EOF
}

while getopts "c:i:dr:t:vhln" option
do
    case $option in
  c) cluster_centers="$OPTARG";;
  i) input="$OPTARG";;
  d) de_novo_cluster=0;;
    r) similarity="$OPTARG";;
  t) threads="$OPTARG";;
    v) verbose=0;;
  l) left_gaps_allowed=0;;
  n) no_overlap=0;;
    h) print_help; exit 0;;
    [?]) print_help; exit 1;;
    esac
done

print_message()
{
    if [ $verbose ]
    then
    echo "`date +%T` $1" >&2
    fi
}

parameters=""
if [ $left_gaps_allowed ]
then
  parameters+=" --left-gaps-allowed "
fi

if [ $no_overlap ]
then
  parameters+=" --no-overlap "
fi

dnaclust_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
print_message "$dnaclust_path"
#exit 1
tempdir=`mktemp -d -p .`
#tempdir="tmpref/"
trap "rm -fr $tempdir" EXIT

#sequences_sorted=`mktemp -p $tempdir`
db_sorted=`basename ${cluster_centers} .fasta`.sorted.fasta
# Reads the sequences from STDIN.
print_message "Reading and sorting the database sequences: $tempdir/${db_sorted}" 
#"$dnaclust_path/fastasort" --random-shuffle > $sequences_sorted 
cat $cluster_centers | "$dnaclust_path/fastasort" > $tempdir/${db_sorted}

print_message "Recruiting from the sequences, using database."
"$dnaclust_path/dnaclust" $parameters -s $similarity -t $threads --no-k-mer-filter -i $input -p $tempdir/${db_sorted} -r | awk '{if (NF > 1) print $0}' > $input.db.clusters