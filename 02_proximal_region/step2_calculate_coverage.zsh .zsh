#!/bin/zsh

# The genome folder stores the fasta files of all candidate genomes.

for a in ../genome/*fasta;
do 
an=`basename -s .fasta $a`
for b in ../genome/*fasta;
do 
bn=`basename -s .fasta $b`
cat ${an}..${bn}.paf | awk '{if($1==$6&&$12>=20) print $6 "\t" $8 "\t" $9}' | bedtools sort -i /dev/stdin | bedtools merge -i /dev/stdin | bedtools coverage -a ../genome/${an}.fasta.fai.5Mb.bed -b /dev/stdin > ${an}..${bn}.ref.cov &
cat ${an}..${bn}.paf | awk '{if($1==$6&&$12>=20) print $1 "\t" $3 "\t" $4}' | bedtools sort -i /dev/stdin | bedtools merge -i /dev/stdin | bedtools coverage -a ../genome/${bn}.fasta.fai.5Mb.bed -b /dev/stdin > ${an}..${bn}.qry.cov

done
done
