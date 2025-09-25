#!/bin/zsh

# The genome folder stores the fasta files of all candidate genomes.

for a in `ls genome/*fasta`
do
for b in `ls genome/*fasta`
do

echo $a
echo $b
aname=`basename -s .fasta $a `
bname=`basename -s .fasta $b `
echo $aname
echo $bname

ref='/filer-dg/agruppen/dg7/fengj/panh2/chr4H/genome/'$aname'.fasta'
qry='/filer-dg/agruppen/dg7/fengj/panh2/chr4H/genome/'$bname'.fasta'

prefix=$aname'..'$bname

minimap2 -t 5 -2 -K 5G -f 0.005 -x asm10 $ref $qry > paf/\$prefix.paf

cd /filer-dg/agruppen/dg7/fengj/panh2/chr4H
done
done
