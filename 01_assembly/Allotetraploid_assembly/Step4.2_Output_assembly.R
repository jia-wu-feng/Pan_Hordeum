#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')


#Function to combine the Hi-C maps from two subgenomes
combine_hic <- function(hap1, hap2, assembly, species="hordeum_bulbosum"){
 assembly_v2 <- assembly
 hic_map_v1_hap1 <- hap1
 hic_map_v1_hap2 <- hap2
 hic_map_v1_hap1$agp[agp_chr != "chrUn"] -> a1
 hic_map_v1_hap2$agp[agp_chr != "chrUn"] -> a2
 a1[, agp_chr := paste0(agp_chr, "_1")]
 a2[, agp_chr := paste0(agp_chr, "_2")]
 a1[, chr := NULL]
 a2[, chr := NULL]
 chrNames(agp=T, species)[, .(chr, agp_chr)][a1, on="agp_chr"] -> a1
 chrNames(agp=T, species)[, .(chr, agp_chr)][a2, on="agp_chr"] -> a2

 c(a1[scaffold != "gap"]$scaffold, a2[scaffold != "gap"]$scaffold) -> s
 s[duplicated(s)] -> s

 a1[s, on="scaffold", scaffold := paste0(scaffold, "_hap1")]
 a2[s, on="scaffold", scaffold := paste0(scaffold, "_hap2")]
 rbind(a1, a2) -> a
 
 hic_map_v1_hap1$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_1"), length, truechr)] -> l1
 hic_map_v1_hap2$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_2"), length, truechr)] -> l2
 rbind(l1, l2) -> l
 l[, offset := cumsum(c(0, length[1:(.N-1)]))]
 l[, plot_offset := cumsum(c(0, length[1:(.N-1)]+1e8))]
 chrNames(agp=T, species)[l, on="agp_chr"] -> l

 copy(assembly_v2$info) -> ai
 ai[!s, on="scaffold"] -> u
 ai[s, on="scaffold"][, scaffold := paste0(scaffold, "_hap1")] -> i1
 ai[s, on="scaffold"][, scaffold := paste0(scaffold, "_hap2")] -> i2
 rbind(u, i1, i2) -> i

 assembly_v2$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> f
 f[scaffold1 %in% s, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]
 f[scaffold2 %in% s, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]

 list(info=i, fpairs=f) -> assembly_hap
 list(agp=a, chrlen=l) -> hic_map
 list(assembly_hap=assembly_hap, hic_map=hic_map)
}

readRDS('hjubatum_assembly_v2.Rds') -> assembly_v2

#import Hi-C maps after editing, v2 -> v3 subgenome change (run the mash for check the distance and change in the excel table)

read_hic_map(rds="hjubatum_hic_map_v1_hap1.Rds", file="hjubatum_hic_map_v1_hap1_edit1.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v3_hap1
saveRDS(hic_map_v3_hap1, file="hjubatum_hic_map_v3_hap1.Rds")

read_hic_map(rds="hjubatum_hic_map_v1_hap2.Rds", file="hjubatum_hic_map_v1_hap2_edit1.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v3_hap2
saveRDS(hic_map_v3_hap2, file="hjubatum_hic_map_v3_hap2.Rds")

#proceed with Hi-C plots and repeat edit cycle if need be



#create Hi-C plots for Subgenome 1
snuc <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hjubatum_hic_map_v3_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v3_hap1_l

#create Hi-C plots for Subgenome 2
snuc <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hjubatum_hic_map_v3_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v3_hap2_l



#readRDS('hjubatum_assembly_v2.Rds') -> assembly_v2
readRDS('hjubatum_hic_map_v3_hap1.Rds') -> hic_map_v3_hap1
readRDS('hjubatum_hic_map_v3_hap2.Rds') -> hic_map_v3_hap2


fasta <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm.fasta'
sink("hjubatum_pseudomolecules_v3_hap1.log")
compile_psmol(fasta=fasta, output="hjubatum_pseudomolecules_v3_hap1", hic_map=hic_map_v3_hap1, assembly=assembly_v2, cores=30)
sink()

sink("hjubatum_pseudomolecules_v3_hap2.log")
compile_psmol(fasta=fasta, output="hjubatum_pseudomolecules_v3_hap2", hic_map=hic_map_v3_hap2, assembly=assembly_v2, cores=30)
sink()


combine_hic(hic_map_v3_hap1, hic_map_v3_hap2, assembly_v2) -> ch2
ch2$hic_map -> hic_map_v3
ch2$assembly -> assembly_v3_hap

nuc <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v3_hap, hic_map=hic_map_v3, nucfile=nuc)->hic_map_v3_l

bin_hic_step(hic=hic_map_v3_l$links, frags=hic_map_v3_l$frags, binsize=1e6, chrlen=hic_map_v3_l$chrlen, cores=14)->hic_map_v3_l$hic_1Mb

f <- "hjubatum_hic_map_v3_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v3_l, file=f, species="hordeum_bulbosum")



