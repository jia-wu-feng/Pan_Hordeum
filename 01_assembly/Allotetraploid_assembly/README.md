# Pan_Hordeum
Code for Article "A genus-wide pangenome of Hordeum"

## Allotetraploid genome assembly (example accession: H. jubatum BCC 2055)

### Step 1 

#### Genome size evaluation 

````shell
zcat ../hifi_reads/*.fastq.gz | jellyfish count /dev/fd/0 -C -o hjubatum -m 71 -t 20 -s 5G

jellyfish histo -h 3000000 -o hjubatum.histo hjubatum
````

````R
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))

library("findGSE")

findGSE(histo="hjubatum.histo", sizek=71, outdir="hjubatum_k71", exp_hom=50)
````

#### Assembly using HiFiasm

Step1.1_Hifiasm.zsh

Convert gfa to fasta and calculate N50

````shell
gfatools gfa2fa hjubatum.p_ctg.gfa | seqtk rename - 'contig_' > hjubatum_hifiasm.fasta

samtools faidx hjubatum_hifiasm.fasta

/filer-dg/agruppen/seq_shared/mascher/code_repositories/tritexassembly.bitbucket.io/shell/n50 hjubatum_hifiasm.fasta > hjubatum_hifiasm.fasta.n50

````

#### Align Hi-C data to assembled unitigs

Step1.2_Hi-C_map.zsh

#### Use the barley high-confidence gene as the maker to build a guide map

Step1.3_Guide_map.zsh

#### create assembly object 

Step1.4_Create_assembly_object.R


### Step 2 

#### Subgenome phasing

Step2_Subgenome_Phasing.R

### Step 3 

#### Output primary assembly

Step3_Output_primary_assembly.R

### Step 4 

#### Ancestral genome-guided subgenome orientation and final output

Step4.1_subgenome_orientation.zsh

Step4.2_Output_assembly.R

##### In allopolyploids without an available ancestral genome, SubPhaser(https://github.com/zhangrengang/SubPhaser) is employed to resolve the chromosomal origins.

#### Output unassembled sequence

````R
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

readRDS('hjubatum_hic_map_v2_hap1.Rds') -> hic_map_v2_hap1
readRDS('hjubatum_hic_map_v2_hap2.Rds') -> hic_map_v2_hap2

hic_map_v2_hap1$agp[!is.na(chr)]$scaffold->h1
hic_map_v2_hap2$agp[!is.na(chr)]$scaffold->h2

hic_map_v2_hap1$agp[,.(scaffold,scaffold_length)][scaffold!="gap"] -> chrh

chrh[,hap:=0]
chrh[scaffold %in% h1, hap:=hap+10]
chrh[scaffold %in% h2, hap:=hap+2]


write.table(chrh,"hjubatum_contig_all.txt",quote=F,row.names=F,sep="\t")

sum(hic_map_v2_hap1$agp[!is.na(chr)]$scaffold_length)
sum(hic_map_v2_hap2$agp[!is.na(chr)]$scaffold_length)


sum(chrh[hap==0]$scaffold_length)

````


````shell
less ../hjubatum_contig_all.txt | awk '{if($3==0) print $1}' > uncontig

seqkit grep -f uncontig ../hjubatum_pseudomolecules_v2_hap1/*_assembly_v2.fasta > hjubatum_unanchor.fasta

````