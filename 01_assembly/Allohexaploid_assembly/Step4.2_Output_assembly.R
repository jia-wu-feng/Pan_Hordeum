#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')


chrNames <- function(agp=F, species="wheat") {
 if(species == "wheat"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 } else if (species == "barley"){
  data.table(alphachr=apply(expand.grid(1:7, "H", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "rye"){
  data.table(alphachr=apply(expand.grid(1:7, "R", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "lolium"){
  data.table(alphachr=as.character(1:7), chr=1:7)->z
 }
 else if (species == "maize"){
  data.table(alphachr=as.character(1:10), chr=1:10)->z
 }
 else if (species == "sharonensis"){
  data.table(alphachr=apply(expand.grid(1:7, "S", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "oats"){
  data.table(alphachr=sub(" ", "", apply(expand.grid(1:21, "M", stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 else if (species == "oats_new"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "C", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 }
 else if (species == "avena_barbata"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:2), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum_4x"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:4), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:28)->z
 }
 else if (species == "hordeum_hexaploid_6x"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:3), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 if(agp){
  rbind(z, data.table(alphachr="Un", chr=0))[, agp_chr := paste0("chr", alphachr)]->z
 }
 z[]
}

wheatchr <- function(agp=F, species="wheat") {
 if(species == "wheat"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 } else if (species == "barley"){
  data.table(alphachr=apply(expand.grid(1:7, "H", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "rye"){
  data.table(alphachr=apply(expand.grid(1:7, "R", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "lolium"){
  data.table(alphachr=as.character(1:7), chr=1:7)->z
 }
 else if (species == "maize"){
  data.table(alphachr=as.character(1:10), chr=1:10)->z
 }
 else if (species == "sharonensis"){
  data.table(alphachr=apply(expand.grid(1:7, "S", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "oats"){
  data.table(alphachr=sub(" ", "", apply(expand.grid(1:21, "M", stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 else if (species == "oats_new"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "C", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 }
 else if (species == "avena_barbata"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:2), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum_4x"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:4), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:28)->z
 }
 else if (species == "hordeum_hexaploid_6x"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:3), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 if(agp){
  rbind(z, data.table(alphachr="Un", chr=0))[, agp_chr := paste0("chr", alphachr)]->z
 }
 z[]
}



#Function to combine the Hi-C maps from two haplotypes
combine_hic <- function(hh,hap1, hap2, hap3, assembly, species="hordeum_hexaploid_6x"){
 assembly_v2 <- assembly
 
 hic_map_v1_hap1 <- hap1
 hic_map_v1_hap2 <- hap2
 hic_map_v1_hap3 <- hap3
 
 hic_map_v1_hap1$agp[agp_chr != "chrUn"] -> a1
 hic_map_v1_hap2$agp[agp_chr != "chrUn"] -> a2
 hic_map_v1_hap3$agp[agp_chr != "chrUn"] -> a3
 
 a1[, agp_chr := paste0(agp_chr, "_1")]
 a2[, agp_chr := paste0(agp_chr, "_2")]
 a3[, agp_chr := paste0(agp_chr, "_3")]
 
 a1[, chr := NULL]
 a2[, chr := NULL]
 a3[, chr := NULL]
 
 chrNames(agp=T, species)[, .(chr, agp_chr)][a1, on="agp_chr"] -> a1 ###### there for ask  
 chrNames(agp=T, species)[, .(chr, agp_chr)][a2, on="agp_chr"] -> a2
 chrNames(agp=T, species)[, .(chr, agp_chr)][a3, on="agp_chr"] -> a3

 c(a1[scaffold != "gap"]$scaffold, a2[scaffold != "gap"]$scaffold, a3[scaffold != "gap"]$scaffold) -> s
 s[duplicated(s)] -> s

 a1[hh[grepl("1", hap) & hap>=10]$scaffold, on="scaffold", scaffold := paste0(scaffold, "_hap1")]
 a2[hh[grepl("2", hap) & hap>=10]$scaffold, on="scaffold", scaffold := paste0(scaffold, "_hap2")]
 a3[hh[grepl("3", hap) & hap>=10]$scaffold, on="scaffold", scaffold := paste0(scaffold, "_hap3")]
 rbind(a1, a2, a3) -> a
 
 hic_map_v1_hap1$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_1"), length, truechr)] -> l1
 hic_map_v1_hap2$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_2"), length, truechr)] -> l2
 hic_map_v1_hap3$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_3"), length, truechr)] -> l3

 rbind(l1, l2, l3) -> l
 l[, offset := cumsum(c(0, length[1:(.N-1)]))]
 l[, plot_offset := cumsum(c(0, length[1:(.N-1)]+1e8))]
 chrNames(agp=T, species)[l, on="agp_chr"] -> l

 copy(assembly_v2$info) -> ai
 ai[!s, on="scaffold"] -> u
 ai[hh[grepl("1", hap) & hap >= 10]$scaffold, on="scaffold"][, scaffold := paste0(scaffold, "_hap1")] -> i1
 ai[hh[grepl("2", hap) & hap >= 10]$scaffold, on="scaffold"][, scaffold := paste0(scaffold, "_hap2")] -> i2
 ai[hh[grepl("3", hap) & hap >= 10]$scaffold, on="scaffold"][, scaffold := paste0(scaffold, "_hap3")] -> i3

 rbind(u, i1, i2, i3) -> i

 assembly_v2$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> f
 f[scaffold1 %in% hh[hap==1234]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), ifelse(runif(.N) > 0.5, "_hap3","_hap4" )))] 
 f[scaffold1 %in% hh[hap==123]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), "_hap3"))]
 f[scaffold1 %in% hh[hap==134]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap3"), "_hap4"))]
 f[scaffold1 %in% hh[hap==124]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), "_hap4"))]
 f[scaffold1 %in% hh[hap==234]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap2","_hap3"), "_hap4"))]
 f[scaffold1 %in% hh[hap==12]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]
 f[scaffold1 %in% hh[hap==13]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap1", "_hap3"))]
 f[scaffold1 %in% hh[hap==14]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap1", "_hap4"))]
 f[scaffold1 %in% hh[hap==23]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap2", "_hap3"))]
 f[scaffold1 %in% hh[hap==24]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap2", "_hap4"))]
 f[scaffold1 %in% hh[hap==34]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap3", "_hap4"))]
 
 f[scaffold2 %in% hh[hap==1234]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), ifelse(runif(.N) > 0.5, "_hap3","_hap4" )))] 
 f[scaffold2 %in% hh[hap==123]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), "_hap3"))]
 f[scaffold2 %in% hh[hap==134]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap3"), "_hap4"))]
 f[scaffold2 %in% hh[hap==124]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), "_hap4"))]
 f[scaffold2 %in% hh[hap==234]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap2","_hap3"), "_hap4"))]
 f[scaffold2 %in% hh[hap==12]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]
 f[scaffold2 %in% hh[hap==13]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap1", "_hap3"))]
 f[scaffold2 %in% hh[hap==14]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap1", "_hap4"))]
 f[scaffold2 %in% hh[hap==23]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap2", "_hap3"))]
 f[scaffold2 %in% hh[hap==24]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap2", "_hap4"))]
 f[scaffold2 %in% hh[hap==34]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap3", "_hap4"))]

 list(info=i, fpairs=f) -> assembly_hap
 list(agp=a, chrlen=l) -> hic_map
 list(assembly_hap=assembly_hap, hic_map=hic_map)
}


readRDS('hparodii_assembly_v2.Rds') -> assembly_v2


#import Hi-C maps after editing, v1 -> v3, change for subgenome assignment
read_hic_map(rds="hparodii_hifiasm_map_v1_hap1.Rds", file="hparodii_hifiasm_map_v1_hap1_edit1.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v3_hap1
saveRDS(hic_map_v3_hap1, file="hparodii_hifiasm_map_v3_hap1.Rds")

read_hic_map(rds="hparodii_hifiasm_map_v1_hap2.Rds", file="hparodii_hifiasm_map_v1_hap2_edit1.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v3_hap2
saveRDS(hic_map_v3_hap2, file="hparodii_hifiasm_map_v3_hap2.Rds")

read_hic_map(rds="hparodii_hifiasm_map_v1_hap3.Rds", file="hparodii_hifiasm_map_v1_hap3_edit1.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v3_hap3
saveRDS(hic_map_v3_hap3, file="hparodii_hifiasm_map_v3_hap3.Rds")


#create Hi-C plots for subgenome 1
snuc <- '../hparodii_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hparodii_hifiasm_map_v3_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v3_hap1_l

#create Hi-C plots for subgenome 2
snuc <- '../hparodii_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hparodii_hifiasm_map_v3_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v3_hap2_l

#create Hi-C plots for subgenome 3
snuc <- '../hparodii_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hparodii_hifiasm_map_v3_hap3.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v3_hap3_l

readRDS('hparodii_hifiasm_map_v3_hap1.Rds') -> hic_map_v3_hap1
readRDS('hparodii_hifiasm_map_v3_hap2.Rds') -> hic_map_v3_hap2
readRDS('hparodii_hifiasm_map_v3_hap3.Rds') -> hic_map_v3_hap3


fasta <- '/filer-dg/agruppen/dg7/fengj/genome/hparodii/hparodii.fasta'

sink("hparodii_pseudomolecules_v3_hap1.log")
compile_psmol(fasta=fasta, output="hparodii_pseudomolecules_v3_hap1", hic_map=hic_map_v3_hap1, assembly=assembly_v2, cores=30)
sink()

sink("hparodii_pseudomolecules_v3_hap2.log")
compile_psmol(fasta=fasta, output="hparodii_pseudomolecules_v3_hap2", hic_map=hic_map_v3_hap2, assembly=assembly_v2, cores=30)
sink()

sink("hparodii_pseudomolecules_v3_hap3.log")
compile_psmol(fasta=fasta, output="hparodii_pseudomolecules_v3_hap3", hic_map=hic_map_v3_hap3, assembly=assembly_v2, cores=30)
sink()


readRDS(file="hparodii_assembly_v2_Hv_guide+HiClift+hap1.Rds") -> hh

hic_map_v3_hap1$agp[agp_chr!="chrUn"][scaffold!="gap"]$scaffold-> hap1
hic_map_v3_hap2$agp[agp_chr!="chrUn"][scaffold!="gap"]$scaffold-> hap2
hic_map_v3_hap3$agp[agp_chr!="chrUn"][scaffold!="gap"]$scaffold-> hap3

hh$scaffold-> hscaffold

for (i in hscaffold) {
  if(i %in% hap1){
     hh[scaffold==i, hap:=1]
  }
  if(i %in% hap2){
     hh[scaffold==i, hap:=2]
  }
  if(i %in% hap3){
     hh[scaffold==i, hap:=3]
  }
  if(i %in% hap1 & i %in% hap2){
     hh[scaffold==i, hap:=12]
  }
  if(i %in% hap1 & i %in% hap3){
     hh[scaffold==i, hap:=13]
  }
  if(i %in% hap2 & i %in% hap3){
     hh[scaffold==i, hap:=23]
  }
  if (i %in% hap1 & i %in% hap2 & i %in% hap3){
     hh[scaffold==i, hap:=123]
  }
}

combine_hic(hh,hic_map_v3_hap1, hic_map_v3_hap2, hic_map_v3_hap3,assembly_v2) -> ch3
ch3$hic_map -> hic_map_v3
ch3$assembly_hap -> assembly_v3_hap

#create inter-chromosomal plot for combined subgenome1-3 Hi-C haplotypes
nuc <- '../hparodii_MboI_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v3_hap, hic_map=hic_map_v3, nucfile=nuc)->hic_map_v3_l

bin_hic_step(hic=hic_map_v3_l$links, frags=hic_map_v3_l$frags, binsize=1e6, chrlen=hic_map_v3_l$chrlen, cores=14)->hic_map_v3_l$hic_1Mb

f <- "hparodii_map_v3_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v3_l, file=f,  species="hordeum_hexaploid_6x")


