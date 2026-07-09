#!/bin/bash

# Blast fasta of candidate genes onto reference genomes (that will be present in the tree)
module load blast-plus

species_list=$1

while read species ; do

# set the necessary variables
ref_genomes_dir=$2
candidates=$3
out_dir=$4

mkdir -p $out_dir

#change working directory
cd $ref_genomes

for f in *fa ; do

blastn -query $candidates -db $f -out $out_dir/candidates_vs_${f/.fa/.out6} -outfmt 6 -num_threads 20 -evalue 1e-10

done

cd -

done < $species_list
