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

##
# Assign contigs that are not place in the guide map to subgenome using Hi-C
##

#read haplotype separation and guide map table
readRDS(file="hjubatum_assembly_v2_haplotype_separation_v0.Rds") -> pq
readRDS(file="hjubatum_assembly_v2_Hv_guide+HiClift.Rds") -> cc

readRDS('hjubatum_assembly_v2.Rds') -> assembly_v2
readRDS('hjubatum_assembly_v2_cov.Rds') -> cov


#find and tabulate links between unplaced and placed con tigs
assembly_v2$fpairs[, .(scaffold=scaffold1, pos=pos1, link=scaffold2)] -> ff
assembly_v2$info[, .(scaffold, length)][ff, on="scaffold"] -> ff
ff[pos <= 2e6 | length - pos <= 2e6] -> ff
ff[scaffold != link] -> ff
ff[, .N, key=.(scaffold, link)] -> fa
fa[N > 1] -> fa
fa[scaffold %in% setdiff(cc$scaffold, pq$contig)] -> fa
fa[link %in% pq[hap %in% 1:2]$contig] -> fa
cc[, .(scaffold, scaffold_chr=chr, scaffold_pos=pos)][fa, on="scaffold"] -> fa
cc[, .(link=scaffold, link_chr=chr, link_pos=pos)][fa, on="link"] -> fa
pq[, .(link=contig, hap)][fa, on="link"] -> fa
fa[scaffold_chr == link_chr] -> fa
fa[, .(n=sum(N)), key=.(scaffold, chr=scaffold_chr, hap)] -> fv
fv[, p := n/sum(n), by=scaffold]
assembly_v2$info[, .(scaffold, length)][fv, on="scaffold"] -> fv
fv[length >= 0][order(-p)][!duplicated(scaffold)][order(-length)] -> hap_lift0
cov[, .(scaffold, cc)][hap_lift0, on="scaffold"] -> hap_lift0

saveRDS(hap_lift0, file="hjubatum_assembly_v2_hap_lift0.Rds") 

#keep contigs unassigned to one haplotype OR contigs assigned to both haplotyped with double coverage
hap_lift0[p >= 0.7] -> hap_lift

#keep onlt contigs >= 300 kb
hap_lift[length >= 3e5] -> hap_lift

#combine both haplotype tables
rbind(pq[, .(scaffold=contig, hap)],  hap_lift[, .(scaffold, hap)]) -> hh

hh[cc, on="scaffold"] -> hh


###check identity and decide the subgenome type (but it's not work good, because the barley almost have same distance for other 2 genome)

aa[,mean(id), key=.(scaffold)] -> chrgeneid
hh[chrgeneid, on="scaffold"] -> chr_sub

chr_sub[,mean(V1), key=.(hap,chr)] -> chr_sub2

for( i  in morexfai[1:7, chr]){
  if(chr_sub2[chr==i&hap==1]$V1 < chr_sub2[chr==i&hap==2]$V1){
     hh[chr==i&hap==1,hap:=3]
     hh[chr==i&hap==2,hap:=4]
     hh[chr==i&hap==3,hap:=2]
     hh[chr==i&hap==4,hap:=1]
  }
}

saveRDS(hh, file="hjubatum_assembly_v2_Hv_guide+HiClift+hap.Rds") 


##
# Construct first Hi-C map
##

#read Hi-C fragment data
f <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp.bed'
read_fragdata(info=assembly_v2$info, file=f)->frag_data

#read subgenome information
readRDS(file="hjubatum_assembly_v2_Hv_guide+HiClift+hap.Rds") -> hh

#Hi-C map for subgenome 1
hh[hap %in% c(1,3), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
  min_nfrag_scaffold=50, max_cM_dist = 1000,
  binsize=3e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap1

saveRDS(hic_map_v1_hap1, file="hjubatum_hic_map_v1_hap1.Rds")

#Hi-C map for subgenome 2
hh[hap %in% c(2,3), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
  min_nfrag_scaffold=50, max_cM_dist = 1000,
  binsize=3e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap2

saveRDS(hic_map_v1_hap2, file="hjubatum_hic_map_v1_hap2.Rds")

#create Hi-C plots for subgenome 1
snuc <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hjubatum_hic_map_v1_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap1_l

#create Hi-C plots for subgenome 2
snuc <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hjubatum_hic_map_v1_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap2_l

#combine both Hi-C maps
combine_hic(hic_map_v1_hap1, hic_map_v1_hap2, assembly_v2) -> ch
ch$hic_map -> hic_map_v1
ch$assembly -> assembly_v2_hap

#create inter-chromosomal plot for combined hap1+hap2 Hi-C haplotypes
nuc <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2_hap, hic_map=hic_map_v1, nucfile=nuc)->hic_map_v1_l

bin_hic_step(hic=hic_map_v1_l$links, frags=hic_map_v1_l$frags, binsize=1e6, chrlen=hic_map_v1_l$chrlen, cores=14)->hic_map_v1_l$hic_1Mb

f <- "hjubatum_hic_map_v1_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v1_l, file=f, species="hordeum_bulbosum")

#write Hi-C map for editing
write_hic_map(rds="hjubatum_hic_map_v1_hap1.Rds", file="hjubatum_hic_map_v1_hap1.xlsx", species="barley")
write_hic_map(rds="hjubatum_hic_map_v1_hap2.Rds", file="hjubatum_hic_map_v1_hap2.xlsx", species="barley")

###### write fasta

readRDS('hjubatum_assembly_v2.Rds') -> assembly_v2
readRDS('hjubatum_hic_map_v1_hap1.Rds') -> hic_map_v1_hap1

readRDS('hjubatum_hic_map_v1_hap2.Rds') -> hic_map_v1_hap2


fasta <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm.fasta'
sink("hjubatum_pseudomolecules_v1_hap1.log")
compile_psmol(fasta=fasta, output="hjubatum_pseudomolecules_v1_hap1", hic_map=hic_map_v1_hap1, assembly=assembly_v2, cores=30)
sink()

sink("hjubatum_pseudomolecules_v1_hap2.log")
compile_psmol(fasta=fasta, output="hjubatum_pseudomolecules_v1_hap2", hic_map=hic_map_v1_hap2, assembly=assembly_v2, cores=30)
sink()

####run alginment for 2 hapltypes 


####Subgenome1


.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source("/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R")

read_paf("hv_hjubatum_pseudomolecules_v1_hap1.paf.gz", save=F, primary_only=T) -> paf #### A_B qry=B ref=A 
paf[, reference := sub("_1", "", reference)]
paf[query == reference, ] -> paf
paf[mapq >= 30] -> paf
paf[alnlen >= 2000] -> paf 
paf[, chr := sub("chr", "", query)]
paf[, chr := sub("H", "", chr)]
paf[, idx := 1:.N]

pdf("hv_hjubatum_v1_hap1_correlation_plots.pdf") 
par(mar=c(5,5,3,1))
par(cex.main=1)
par(cex.lab=1)
par(cex.axis=1)
lapply(1:7, function(i){
 paf[chr == i] -> pafL
 pafL[chr == i, plot(las=1, bty='l',type='n', 0, xlab="hjubatum hap1 CCS (Mb)", ylab="FB19 Hap1 CCS (Mb)", xlim=c(0,max(query_end))/1e6, ylim=c(0,max(reference_end))/1e6, main=paste0(i, "H"))] ### xlab = qry = B and ylab = ref = A  
 pafL[orientation == 1, lines(c(query_start/1e6, query_end/1e6), c(reference_start/1e6, reference_end/1e6), col="#000000", lwd=2), by=idx]
 pafL[orientation == -1, lines(c(query_start/1e6, query_end/1e6), c(reference_end/1e6, reference_start/1e6), col="#000000", lwd=2), by=idx]
})
dev.off()



####Subgenome2


.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source("/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R")

read_paf("hv_hjubatum_pseudomolecules_v1_hap2.paf.gz", save=F, primary_only=T) -> paf #### A_B qry=B ref=A 
paf[, reference := sub("_1", "", reference)]
paf[query == reference, ] -> paf
paf[mapq >= 30] -> paf
paf[alnlen >= 2000] -> paf 
paf[, chr := sub("chr", "", query)]
paf[, chr := sub("H", "", chr)]
paf[, idx := 1:.N]

pdf("hv_hjubatum_v1_hap2_correlation_plots.pdf") 
par(mar=c(5,5,3,1))
par(cex.main=1)
par(cex.lab=1)
par(cex.axis=1)
lapply(1:7, function(i){
 paf[chr == i] -> pafL
 pafL[chr == i, plot(las=1, bty='l',type='n', 0, xlab="hjubatum hap2 CCS (Mb)", ylab="FB19 Hap1 CCS (Mb)", xlim=c(0,max(query_end))/1e6, ylim=c(0,max(reference_end))/1e6, main=paste0(i, "H"))] ### xlab = qry = B and ylab = ref = A  
 pafL[orientation == 1, lines(c(query_start/1e6, query_end/1e6), c(reference_start/1e6, reference_end/1e6), col="#000000", lwd=2), by=idx]
 pafL[orientation == -1, lines(c(query_start/1e6, query_end/1e6), c(reference_end/1e6, reference_start/1e6), col="#000000", lwd=2), by=idx]
})
dev.off()


#import Hi-C maps after editing, v1 -> v2
read_hic_map(rds="hjubatum_hic_map_v1_hap1.Rds", file="hjubatum_hic_map_v1_hap1_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap1
saveRDS(hic_map_v2_hap1, file="hjubatum_hic_map_v2_hap1.Rds")

read_hic_map(rds="hjubatum_hic_map_v1_hap2.Rds", file="hjubatum_hic_map_v1_hap2_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap2
saveRDS(hic_map_v2_hap2, file="hjubatum_hic_map_v2_hap2.Rds")

#proceed with Hi-C plots and repeat edit cycle if need be



#create Hi-C plots for subgenome 1
snuc <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hjubatum_hic_map_v2_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v2_hap1_l

#create Hi-C plots for subgenome 2
snuc <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="hjubatum_hic_map_v2_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v2_hap2_l



#readRDS('hjubatum_assembly_v2.Rds') -> assembly_v2
readRDS('hjubatum_hic_map_v2_hap1.Rds') -> hic_map_v2_hap1
readRDS('hjubatum_hic_map_v2_hap2.Rds') -> hic_map_v2_hap2


fasta <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm.fasta'
sink("hjubatum_pseudomolecules_v2_hap1.log")
compile_psmol(fasta=fasta, output="hjubatum_pseudomolecules_v2_hap1", hic_map=hic_map_v2_hap1, assembly=assembly_v2, cores=30)
sink()

sink("hjubatum_pseudomolecules_v2_hap2.log")
compile_psmol(fasta=fasta, output="hjubatum_pseudomolecules_v2_hap2", hic_map=hic_map_v2_hap2, assembly=assembly_v2, cores=30)
sink()


combine_hic(hic_map_v2_hap1, hic_map_v2_hap2, assembly_v2) -> ch2
ch2$hic_map -> hic_map_v2
ch2$assembly -> assembly_v2_hap

nuc <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm_MboI_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2_hap, hic_map=hic_map_v2, nucfile=nuc)->hic_map_v2_l

bin_hic_step(hic=hic_map_v2_l$links, frags=hic_map_v2_l$frags, binsize=1e6, chrlen=hic_map_v2_l$chrlen, cores=14)->hic_map_v2_l$hic_1Mb

f <- "hjubatum_hic_map_v2_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v2_l, file=f, species="hordeum_bulbosum")




