#!/bin/zsh

### The folder(/filer-dg/agruppen/dg7/fengj/genome/hparodii/hic/split_subgenome) contains the known ancestral genome (Hordeum flexuosum, Hordeum roshevitzii and Hordeum pubiflorum) and the preliminary V1 assembly.

for i in chr1H chr2H chr3H chr4H chr5H chr6H chr7H;
do
echo $i
name=`basename $i`
echo $name
mkdir $name
cd $name

echo $name > 1

date

for n in `ls /filer-dg/agruppen/dg7/fengj/genome/hparodii/hic/split_subgenome/genome/*`;do echo $n;cname=`basename $n` ;/opt/Bio/seqkit/0.9.1/bin/seqkit grep -f 1 $n > $cname.1.fa; done

for m in `ls /filer-dg/agruppen/dg7/fengj/genome/hparodii/hic/split_subgenome/genome/*`;do echo $m;qname=`basename -s .fasta $m`; echo $qname; sed -i '1c >'$qname'' $qname.fasta.1.fa; done

cat *.fasta.1.fa > ${name}.fa 

bgzip -@ 20 ${name}.fa 

samtools faidx ${name}.fa.gz

/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/pangraph/mash/mash-Linux64-v2.3/mash triangle ${name}.fa.gz > ${name}.stat.txt

cd /filer-dg/agruppen/dg7/fengj/genome/hparodii/hic/split_subgenome

done
