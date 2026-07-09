#!/bin/bash

# Generate consensus sequence for each orthogroup and use that to blast references (that will be present on the tree)
module load emboss
module load mafft
module load blast-plus

OUTDIR_MAFFT=$1
REFDIR=$2
num_threads=$3

for f in $OUTDIR_MAFFT/*full.fa ; do

count=$(grep -c ">" $f)

if [ $count -eq 1 ] ; then 
cp $f ${f/.fa/_cons.fa}
echo "$f - single sequence"

else

mafft --thread $num_threads ${f} > ${f/.fa/_mafft.fa}
cons -sequence ${f/.fa/_mafft.fa} -outseq ${f/.fa/_mafft_cons.fa}

fi

done


for ref in $REFDIR/*fa ; do

ref_name=${ref/.fa/}

makeblastdb -in $ref -dbtype nucl 

for f in *cons.fa ; do

blastn -query ${f} -db $ref -out ${f/.fa}_vs_${ref_name}.out6 -outfmt 6 -num_threads $num_threads -evalue 1e-10

done
