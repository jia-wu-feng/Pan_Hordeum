# Pan_Hordeum
Code for Article "A genus-wide pangenome of Hordeum"

## Allohexaploid genome assembly (example accession: H. hparodii BCC2025)

### Step 1 

#### Genome size evaluation 

````shell
zcat ../hifi_reads/*.fastq.gz | jellyfish count /dev/fd/0 -C -o hjubatum -m 71 -t 20 -s 5G

jellyfish histo -h 3000000 -o hparodii.histo hparodii
````

````R
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))

library("findGSE")

findGSE(histo="hparodii.histo", sizek=71, outdir="hparodii_k71", exp_hom=50)
````

#### Assembly using HiFiasm

Step1.1_Hifiasm.zsh

Convert gfa to fasta and calculate N50

````shell
gfatools gfa2fa hparodii.p_ctg.gfa | seqtk rename - 'contig_' > hparodii.fasta

samtools faidx hparodii.fasta

/filer-dg/agruppen/seq_shared/mascher/code_repositories/tritexassembly.bitbucket.io/shell/n50 hparodii.fasta > hparodii.fasta.n50

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

##### In allopolyploids without an available ancestral genome, SubPhaser(https://github.com/zhangrengang/SubPhaser) is employed to split the subgenome chromosomal.

#### Output unassembled sequence

````R
readRDS('hparodii_hifiasm_map_v3_hap1.Rds') -> hic_map_v3_hap1
readRDS('hparodii_hifiasm_map_v3_hap2.Rds') -> hic_map_v3_hap2
readRDS('hparodii_hifiasm_map_v3_hap3.Rds') -> hic_map_v3_hap3


hic_map_v3_hap1$agp[!is.na(chr)]$scaffold->h1
hic_map_v3_hap2$agp[!is.na(chr)]$scaffold->h2
hic_map_v3_hap3$agp[!is.na(chr)]$scaffold->h3


hic_map_v3_hap1$agp[,.(scaffold,scaffold_length)][scaffold!="gap"] -> chrh

chrh[,hap:=0]
chrh[scaffold %in% h1, hap:=hap+100]
chrh[scaffold %in% h2, hap:=hap+20]
chrh[scaffold %in% h3, hap:=hap+3]


write.table(chrh,"hparodii_contig_all.txt",quote=F,row.names=F,sep="\t")

sum(hic_map_v3_hap1$agp[!is.na(chr)]$scaffold_length)
sum(hic_map_v3_hap2$agp[!is.na(chr)]$scaffold_length)
sum(hic_map_v3_hap3$agp[!is.na(chr)]$scaffold_length)

sum(chrh[hap==0]$scaffold_length)

````


````shell
less ../hparodii_contig_all.txt | awk '{if($3==0) print $1}' > uncontig

seqkit grep -f uncontig ../hparodii_pseudomolecules_v2_hap1/*_assembly_v2.fasta > hparodii_unanchor.fasta

````