# Pan_Hordeum
Code for Article "A genus-wide pangenome of Hordeum"

## Identification of SVs

### Step 1 SV identification based on wfmash and SyRI

shell
```
for i in genome/*.fasta;
do
echo $i
name=`basename $i`
echo $name
mkdir $name
cd $name

ref='/filer-dg/agruppen/dg7/fengj/genome/syri_wfmash_f/reference/MorexV3_pseudomolecules.fasta'
qry='/filer-dg/agruppen/dg7/fengj/genome/syri_wfmash_f/genome/$name'

prefix='MorexV3_$name'

wfmash -p 90 -l 0 -t 10 $ref $qry > $prefix.paf 2> $prefix.err 

syri -c ${prefix}.paf -F P -r $ref -q $qry --prefix $prefix -f --nc 7 --nosnp 

cd /filer-dg/agruppen/dg7/fengj/genome/syri_wfmash_f
done
```

### Step 2 Extraction of large inversions (>2 Mb)

shell
```
for i in *.vcf; 
do 
echo $i; 
cat $i | grep "Parent=.;" | grep "<INV>" | awk '{print $1 ";" $2 ";" $3 ";" $8}' | awk -F ';' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 }' | sed s/END=//g | sed s/ChrB=//g | sed s/StartB=//g | sed s/EndB=//g | awk '{if($4-$2>=2000000) print $0}' > $i.INV.txt; 
done

```

### Step 3. Identification of non-proximal large inversion maps

shell
```
for i in *INV.txt; do echo $i; name=`basename -s .INV.bed $i`; echo $name; cat $i | awk '{print $1 "\t" $2 "\t" $4 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" "'${name}'"}' >> all.INV.bed; done

bedtools coverage -a all.INV.bed -b Hordeum_vulgare_proximal_region.bed | awk '{if($12<=0.8) print $0}' > large_INV_map.bed

```
