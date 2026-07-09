# Pan_Hordeum
Code for Article "A genus-wide pangenome of Hordeum"

## Detection of introgressed regions using gene tree discordance

Introgressed regions were detected by comparing local sliding-window gene trees with the genome-wide reference phylogenetic tree. Single-copy orthologs identified by OrthoFinder were used for both reference tree construction and local tree reconstruction.

### Step 1. Extraction of single-copy orthologous CDS sequences

Single-copy orthologous gene groups were obtained from the OrthoFinder results. Because OrthoFinder was performed using protein sequences, the corresponding CDS sequences were extracted from genome annotation files according to the single-copy gene IDs. Each orthologous gene group was stored in an independent folder, and each FASTA file was named using the format `OG_name___gene_name`.

```
for g in ``;
do
 echo $g
 gname=`basename $g `
 echo $gname
 cd $gname
 cat ${gname}*.fa > all_${gname}.fasta
 mafft --auto all_${gname}.fasta > all_${gname}.al
done

cat */all_*al > all_alignment.fasta

seqkit split -i all_alignment.fasta
```

### Step 2. Construction of the genome-wide reference tree

Aligned single-copy genes were grouped by species according to `species.txt`. The aligned CDS sequences were concatenated following the order of orthologous groups listed in `OG_name.txt`, generating one concatenated sequence for each species.

```
for name in `cat species.txt`;
do
echo $i
name=`basename $i .fa`
echo $name

for g in `cat OG_name.txt`;
do
 echo $g
 gname=`basename $g `
 echo $gname
 cat ${name}/${gname}*.fa >> $name.singlegene.cds.fa
done
done

for i in `ls *.singlegene.cds.fa`; 
do 
echo $i
name=`basename -s .singlegene.cds.fa $i`
echo $name
cat $name.singlegene.cds.fa | grep -v ">" > ${name}_comb.fasta
sed -i ':a;N;$!ba;s/\n//g' ${name}_comb.fasta
sed -i '1i >'$name'' ${name}_comb.fasta
done

cat *_comb.fasta > species_alignment.fasta
```

### Step 3. Filtering of poorly aligned regions and reference tree inference

Poorly aligned regions were removed using Gblocks. The filtered concatenated alignment was then used to construct the genome-wide reference tree with IQ-TREE. `Sorghum_bicolor` was used as the outgroup.

```
Gblocks species_alignment.fasta -t=extractedp -b4=5 -b5=h

iqtree2 -s species_alignment.fasta-gb -m MFP -B 1000 -o Sorghum_bicolor -T AUTO
```

### Step 4. Generation of sliding windows based on gene order in Hordeum vulgare Morex

Single-copy orthologous genes were ordered according to their genomic positions in `Hordeum_vulgare_Morex`. Sliding windows were generated using 25 single-copy genes per window with a step size of 2 genes.

```
for i in chr1H chr2H chr3H chr4H chr5H chr6H chr7H; 
do 
echo $i
less Hordeum_vulgare_Morex_single_copy_gene_list.bed | grep $i | awk '{print $0}' > list.$i
done

for chr in chr1H chr2H chr3H chr4H chr5H chr6H chr7H; 
do 
num=`wc -l list.$chr | awk '{print $1}'`
num1=`awk '{ print int(('$num'-25)/2)}' 1`
echo $num1

for i in {0..$num1}; 
do 
mkdir $chr.$i
n5i=`expr $i \* 2 + 25`
echo $n5i
cd $chr.$i 

head -n $n5i list.$chr | tail -n 25 | awk '{print $5}' > f_50.name

cd /filer-dg/agruppen/dg7/fengj/panh2/single_copy_tree2/tree_methods_w25_s2
done
done
```

### Step 5. Extraction and concatenation of sequences for each sliding window

For each sliding window, CDS sequences of the 25 single-copy genes were extracted for all species and concatenated into one sequence per species.

```
for chr in chr1H chr2H chr3H chr4H chr5H chr6H chr7H; 
do 
num=`wc -l list.$chr | awk '{print $1}'`
num1=`awk '{ print int(('$num'-25)/2)}' 1`
echo $num1

for i in {0..$num1}; 
do 
n5i=`expr $i \* 2 + 25`
echo $n5i
cd $chr.$i 

for i in `cat species.txt`; 
do 
echo $i
name=`basename -s .fa $i`
echo $name

for g in `cat f_50.name`;
do 
echo $g
cat /filer-dg/agruppen/dg7/fengj/panh2/single_copy_tree2/single_copy_tree/$name/${g}___*.fa >> $i.single_copy.fasta
done
done

for i in `cat species.txt`; 
do 
echo $i
name=`basename -s .fa $i`
echo $name
cat $name.fa.single_copy.fasta | grep -v ">" > ${name}_comb.fasta
sed -i ':a;N;$!ba;s/\n//g' ${name}_comb.fasta
sed -i '1i >'$name'' ${name}_comb.fasta
done

cd /filer-dg/agruppen/dg7/fengj/panh2/single_copy_tree2/tree_methods_w25_s2
done
done
```

### Step 6. Construction of local phylogenetic trees for each window

For each sliding window, concatenated sequences from all species were aligned using MAFFT. Conserved alignment blocks were retained using Gblocks, and local phylogenetic trees were inferred using FastTree.

```
for chr in chr1H chr2H chr3H chr4H chr5H chr6H chr7H; 
do 
num=`wc -l tree_list/list.$chr | awk '{print $1}'`
num1=`awk '{ print int(('$num'-25)/2)}' 1`
echo $num1

for i in {0..$num1}; 
do 
n5i=`expr $i \* 2 + 25`
echo $n5i
cd $chr.$i 

cat *_comb.fasta > 18_out_chr1
mafft --thread 30 --auto 18_out_chr1.fasta > 18_out_chr1 
Gblocks 18_out_chr1 -t=extractedp -b4=5 -b5=h
FastTree 18_out_chr1-gb > $chr.$i.tree 

ls $chr.$i.tree

cd /filer-dg/agruppen/dg7/fengj/panh2/single_copy_tree2/tree_methods_w25_s2
done
done
```

### Step 7. Rooting local trees with Sorghum bicolor

All local sliding-window trees were rooted using `Sorghum_bicolor` as the outgroup.

```
for i in *tree
do
name=`basename ${i} .tree`
ete3 mod --outgroup Sorghum_bicolor -t ${i} > ${name}.rt &
sleep 0.1s
done
```

### Step 8. Assignment of genomic coordinates to each sliding window

The genomic range covered by each sliding window was extracted according to gene positions in `Hordeum_vulgare_Morex`.

```
for chr in chr1H chr2H chr3H chr4H chr5H chr6H chr7H; 
do 
num=`wc -l list.$chr | awk '{print $1}'`
num1=`awk '{ print int(('$num'-25)/2)}' 1`
echo $num1

for i in {0..$num1}; 
do 
n5i=`expr $i \* 2 + 25`
echo $n5i

head -n $n5i ../list.$chr | tail -n 25 | awk '{print $0}' > $chr.$i.bed
bedtools merge -d 100000000000000000 -i $chr.$i.bed > $chr.$i.merge.bed

cat $chr.$i.merge.bed | awk '{print "'$chr.$i'.rt" "\t"  $1 "\t" $2 "\t" $3}' >> list.bed
done
done

sort list.bed > sort.list.bed
```

### Step 9. Comparison between local trees and the reference tree

The topological discordance between each local tree and the genome-wide reference tree was quantified using ETE3. The normalized Robinson–Foulds distance was calculated with a minimum branch support threshold of 90.

```
ete3 compare -r all_hordeum.fasta-gb.treefile -t *.rt --min_support_src 90 > tree.scores

sort tree.scores > tree.sort.scores
```

### Step 10. Integration of tree discordance scores with genomic positions

The genomic coordinates of each sliding window were joined with the corresponding tree comparison scores.

```
join sort.list.bed tree.sort.scores -1 1 -2 1 | \
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 }' \
> sort_table_sorce
```

### Step 11. Identification of candidate introgressed regions

Sliding windows with high tree discordance were filtered using the nRF threshold. Adjacent discordant windows within 10 Mb were merged, and regions supported by more than five overlapping discordant windows were defined as candidate introgressed regions.

```
less scoreabove0.21.bed | \
awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' | \
bedtools sort -i /dev/stdin | \
bedtools merge -d 10000000 -i /dev/stdin -c 4 -o count,collapse | \
awk '{if($4>5) print $0}' \
> potatial_events.bed
```

### Step 12. Inference of introgression origin

For each candidate introgressed region, local phylogenetic trees were reconstructed using the same procedure as the genome-wide reference tree. The topology of each local tree was compared with the reference tree to infer the potential donor lineage and determine the origin of large introgressed fragments.