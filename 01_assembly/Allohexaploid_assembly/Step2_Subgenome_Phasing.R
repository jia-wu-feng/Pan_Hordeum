#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')
#read assembly object. This is the hifiasm unitig assembly after breaking chimeras based on coverage drops
readRDS('hparodii_assembly_v2.Rds') -> assembly_v2

##
# Positions contigs based on guide map
##

#read alignment of MorexV3 HC genes 
fread('./hparodii_MorexV3_HC_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))

#read positions of MorexV3 genes
fread('/filer-dg/agruppen/seq_shared/mascher/morex_v3_annotation_200818/Hv_Morex.pgsb.Jul2020.HC_mrna_pos.txt')->p
setnames(p, c("chr", "chr_start", "chr_end", "transcript"))

#keep only genes aligned for >= 90 % of their length with 90 % identity, allow up to three alignment (one per subgenome), modify for allohexaploid
p[x, on="transcript"] -> px
px[alnlen >= 90 & id >= 90]->a
px[a[, .N, key=transcript][N <= 3]$transcript, on="transcript"]->aa

#convert original contig positions to corrected (assembly_v2) positions 
assembly_v2$info[, .(contig=orig_scaffold, start=orig_start, orig_start, scaffold)][aa, on=c("contig", "start"), roll=T] -> aa

#chromosome assignment
aa[, .N, key=.(chr, scaffold)][, p := N/sum(N), by=scaffold][order(-p)][!duplicated(scaffold)] -> cc

#keep contigs with at least 4 aligned genes, 75 % of aligned are from the major chromosome
cc[N >= 4 & p >= 0.75] -> cc

#check how much of the assembly can be assigned to chromosomes
assembly_v2$info[, .(scaffold, scaffold_length=length)][cc, on="scaffold"] -> cc
cc[scaffold_length>=1e6]->cc 
sum(cc$scaffold_length)/1e6

####10066.21 


#get the approximate chromosome positions (median of alignment coordinates)
aa[cc[, .(scaffold, chr)], on=c("chr", "scaffold")][, .(pos=median(as.numeric(chr_start)), pos_mad=mad(as.numeric(chr_start))), key=scaffold][cc, on="scaffold"] -> cc
cc[, mr := pos_mad / scaffold_length]  
setorder(cc, chr, pos) 
cc[, chr_idx := 1:.N, by=chr]
cc[, agp_pos := c(0, cumsum(scaffold_length[-.N])), by=chr]
#save results
saveRDS(cc, file="hparodii_assembly_v2_Hv_guide.Rds")

##
# Positioning additional contigs by Hi-C
## 

#load guide map positions
readRDS('hparodii_assembly_v2_Hv_guide.Rds') -> cc
#get Hi-C links in the terminal 2 Mb of each scaffold
assembly_v2$fpairs[, .(scaffold=scaffold1, pos=pos1, link=scaffold2)] -> ff
assembly_v2$info[, .(scaffold, length)][ff, on="scaffold"] -> ff
ff[pos <= 2e6 | length - pos <= 2e6] -> ff
ff[scaffold != link] -> ff
ff[, .N, key=.(scaffold, link)] -> fa
#exclude scaffold-link paris with only a single Hi-C pair
fa[N > 1] -> fa
#add guide map positions
cc[, .(scaffold, scaffold_chr=chr, scaffold_pos=pos)][fa, on="scaffold"] -> fa
cc[, .(link=scaffold, link_chr=chr, link_pos=pos)][fa, on="link"] -> fa
#get approximate position based on Hi-C links
fa[!is.na(link_chr), .(n=sum(N), pos=weighted.mean(link_pos, N)), key=.(scaffold, scaffold_chr, scaffold_pos, link_chr)] -> fv
fv[, p := n/sum(n), by=scaffold]
assembly_v2$info[, .(scaffold, length)][fv, on="scaffold"] -> fv
#get chromosome assignment for scaffolds without position in guide map
fv[is.na(scaffold_chr)][order(-p)][!duplicated(scaffold)][order(-length)] -> cc_lift0

#exclude contigs shorter than 300 kb
cc_lift0[length >= 3e5] -> cc_lift
#merge guide map and Hi-C lift tables; write output
rbind(cc[, .(scaffold, chr, pos)], cc_lift[, .(scaffold, chr=link_chr, pos)]) -> a
a[assembly_v2$info[, .(scaffold, length)], on="scaffold"] -> a
saveRDS(a, file="hparodii_assembly_v2_Hv_guide+HiClift.Rds")

#check how much of the assembly is placed
 sum(a[!is.na(chr)]$length)
#[1] 10807155785
 sum(a$length)
#[1] 10897711776
 10807155785/10897711776
#[1] 0.9916904

##
# subgenome phasing
##

#read lengths of Morex chromosomes
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/200416_MorexV3_pseudomolecules.fasta.fai', sel=1:2, col.names=c("chr", "len"))->morexfai

#read centromere position in Morex 
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/MorexV3_centromere_positions.tsv')->cen
setnames(cen, c("chr", "cen_pos"))

#read scaffold posiitions and add them to the Hi-C link table
readRDS(file="hparodii_assembly_v2_Hv_guide.Rds") -> cc
assembly_v2$fpairs[, .N, key=.(ctg1=scaffold1, ctg2=scaffold2)] -> v
cc[, .(ctg1=scaffold, chr1=chr, pos1=pos)][v, on="ctg1"] -> v
cc[, .(ctg2=scaffold, chr2=chr, pos2=pos)][v, on="ctg2"] -> v

#run a PCA on the intra-chromosomal matrices
rbindlist(mclapply(mc.cores=7, morexfai[1:7, chr], function(j){
 #fill in empty scaffold pairs with 0
 setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
 dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
 melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
 x[, l := log10(0.1 + N)]
 cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
 cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
 morexfai[x, on="chr"] -> x
 cen[x, on="chr"] -> x
 x[, end_dist1 := pmin(pos1, len - pos1)]
 x[, end_dist2 := pmin(pos2, len - pos2)]
 x[, cen_dist1 := abs(cen_pos - pos1)]
 x[, cen_dist2 := abs(cen_pos - pos2)]
 x[, ldist := abs(pos1 - pos2)]
 #get "normalized" Hi-C ounts by removing the factor linear distance between loci, distance from centromere and distance from chromosome end (all log scaled)
 x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]

 #convert to matrix
 dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
 y[, ctg1 := NULL]
 #run PCA on correlation matrix
 prcomp(cor(y), scale=T, center=T)->pca
 #get first four eigenvector
 data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
 setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp
 pp[, chr := j]
})) -> pp

assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp, on="contig"] -> pp
pp[, idx := 1:.N]
#save results
saveRDS(pp, file="hparodii_assembly_v2_haplotype_separation_HiC.Rds")

#read coverage information and convert to correct assembly coordinates
fread(cmd="grep '^S' ../hparodii.p_ctg.noseq.gfa | cut -f 2,5 | tr ':' '\\t' | cut -f 1,4", head=F) -> cov
setnames(cov, c("contig", "cc"))
cov[, contig := sub('l$', '', sub("ptg0*", "contig_", contig))]
assembly_v2$info[, .(contig=orig_scaffold, scaffold)][cov, on=c("contig")] -> cov
setorder(cc[, .(scaffold, scaffold_length, chr, pos)][cov, on="scaffold"], chr, pos) -> cov
cov[, idx := 1:.N]
saveRDS(file="hparodii_assembly_v2_cov.Rds", cov)


cov[, .(contig=scaffold, cc)][pp, on="contig"] -> pp

pp[cc<45, col := "red"]
pp[cc>=45&cc<80,col:="blue"]
pp[cc>=80&cc<120,col:="black"]
pp[cc>=120,col:="orange"]


pdf("hparodii_haplotype_separation_HiC_all_pc1_2.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pp[chr == i, plot(PC2, PC1, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
  legend("topleft", pch=19, bty='n', legend=c("1X", "2X","3X","4X"), col=c("red", "blue","black","orange"))
 })

dev.off()


pdf("hparodii_haplotype_separation_HiC_all.pdf", height=4, width=10)
par(mfrow=c(2,4))
lapply(c("PC1", "PC2", "PC3", "PC4","cc"), function(j){
 lapply(morexfai[1:7, chr], function(i){
  pp[chr == i, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chr[1]), ylab=j,
            xlab="Hv syntenic position [Mb]")]
  pp[chr == i, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(get(j), get(j))), by=idx]
 })
 plot(axes=F, xlab="", ylab="", 0, type='n')
})
dev.off()

saveRDS(pp, file="hparodii_assembly_v2_haplotype_separation_HiC_all.Rds")


readRDS(file="hparodii_assembly_v2_haplotype_separation_HiC_all.Rds") -> pp4

### Classifying PCA results based on k-means

pp4[,cluster:=0]
primer_cluster <- function(pp4, chrj, center){
data.frame(pp4[chr==chrj,.(PC1,PC2)])-> mm
kmeans(mm,centers=center, iter.max = 30, nstart = 30)$cluster->q1
data.table(pp4[chr==chrj])-> mm1
mm1[,cluster:=q1]
rbind(pp4[chr!=chrj], mm1)->pp4
return(pp4) 
}

primer_cluster(pp4,"chr1H",3)-> pp4
primer_cluster(pp4,"chr2H",3)-> pp4
primer_cluster(pp4,"chr3H",3)-> pp4
primer_cluster(pp4,"chr4H",3)-> pp4
primer_cluster(pp4,"chr5H",3)-> pp4
primer_cluster(pp4,"chr6H",3)-> pp4
primer_cluster(pp4,"chr7H",3)-> pp4

pp[,cluster:=0]->pp1
rbind( pp1[!pp4$contig, on="contig"][, PC1 := 0 ][, PC2 := 0 ][, PC3 := 0 ][, PC4 := 0 ], pp4)->pq

#manually define cuts between haplotype for each chromosome
#2x contigs get "haplotype 3", i.e. present in both 1 and 2
pq[cc>=45 & cc<80, cluster := 22]
pq[cc>=80 & cc<120, cluster := 33]
pq[cc>=120, cluster :=44]

pq[, col := "black"]
pq[cluster == 1, col := "red"]
pq[cluster == 2, col := "blue"]
pq[cluster == 3, col := "purple"]

pq[cluster == 22 |cluster == 33 | cluster == 44, col :="black"]


setorder(pq, chr, pos)->pq
pq[, idx := 1:.N]

#########important
saveRDS(pq, file="hparodii_assembly_v2_manual.Rds")



pdf("hparodii_haplotype_separation_HiC_PC1_2_cluster.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pq[chr == i, plot(PC2, PC1, col=cluster, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
  legend("topleft", pch=19, bty='n', legend=c("1", "2","3","4"), col=c(1,2,3,4))
 })

dev.off()

### Calculate chromosome size after phasing

cal_chrlen <- function(pq, chrj,cov1){
  hap1=sum(pq[chr == chrj & cluster==1 & cc<=cov1]$contig_length)+sum(pq[chr == chrj & cluster==1 & cc>=cov1]$contig_length)*2
  hap2=sum(pq[chr == chrj & cluster==2 & cc<=cov1]$contig_length)+sum(pq[chr == chrj & cluster==2 & cc>=cov1]$contig_length)*2
  hap3=sum(pq[chr == chrj & cluster==3 & cc<=cov1]$contig_length)+sum(pq[chr == chrj & cluster==3 & cc>=cov1]$contig_length)*2
  hap4=sum(pq[chr == chrj & cluster==4 & cc<=cov1]$contig_length)+sum(pq[chr == chrj & cluster==4 & cc>=cov1]$contig_length)*2
  print(hap1)
  print(hap2)
  print(hap3)
  print(hap4)
}

cal_chrlen(pq,"chr1H",45)
cal_chrlen(pq,"chr2H",45)
cal_chrlen(pq,"chr3H",45)
cal_chrlen(pq,"chr4H",45)
cal_chrlen(pq,"chr5H",45)
cal_chrlen(pq,"chr6H",45)
cal_chrlen(pq,"chr7H",45)



###cluster Re- PCA 

Repca <- function(pq, chrj, cov,v, cl, plevel=0)
{
pq=pq
j=chrj
setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
x[ctg1 %in%  pq[cluster==cl ]$contig & ctg2 %in% pq[cluster==cl ]$contig] -> x
dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
x[, l := log10(0.1 + N)]
cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
morexfai[x, on="chr"] -> x
cen[x, on="chr"] -> x
x[, end_dist1 := pmin(pos1, len - pos1)]
x[, end_dist2 := pmin(pos2, len - pos2)]
x[, cen_dist1 := abs(cen_pos - pos1)]
x[, cen_dist2 := abs(cen_pos - pos2)]
x[, ldist := abs(pos1 - pos2)]
x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]
dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
y[, ctg1 := NULL]
prcomp(cor(y), scale=T, center=T)->pca
data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp1c1
pp1c1[, chr := j]
assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp1c1, on="contig"] -> pp1c1
pp1c1[, idx := 1:.N]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[cc<32, col := "red"]
pp1c1[cc>=32&cc<64,col:="blue"]
pp1c1[cc>=64&cc<96,col:="black"]
pp1c1[cc>=96,col:="orange"]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[PC1>plevel, cluster := cl]
pp1c1[PC1<plevel, cluster := 4]
pdf(paste("haplotype_separation_HiC_pp",chrj,"c",cl,".pdf",sep=""), height=8, width=8)
plot(pp1c1$PC2, pp1c1$PC1, col=pp1c1$col, pch=20, las=1,bty='l', type="p", main=paste(chrj," cluster",cl,sep=""), ylab="PC1",xlab="PC2")
dev.off()
rbind(pq[!pp1c1$contig, on="contig"], pp1c1)-> pq1
print(sum(pq1[chr == chrj & cluster==1]$contig_length))
print(sum(pq1[chr == chrj & cluster==2]$contig_length))
print(sum(pq1[chr == chrj & cluster==3]$contig_length))
print(sum(pq1[chr == chrj & cluster==4]$contig_length))
pp1c1[, idx := 1:.N]
pdf(paste("haplotype_separation_HiC_pp",chrj,"c",cl,"_repca.pdf",sep=""), height=4, width=10)
lapply(c("PC1", "PC2", "PC3", "PC4"), function(j){
 pp1c1[chr==chrj, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chrj), ylab=j,
            xlab="Hv syntenic position [Mb]")]
 pp1c1[chr==chrj, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(get(j), get(j))), by=idx]
})
dev.off()
return(pq1) 
}

Repcafor2 <- function(pq, chrj, cov,v, cl,cl2, plevel=0)
{
pq=pq
j=chrj
setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
x[ctg1 %in%  pq[cluster==cl | cluster==cl2 ]$contig & ctg2 %in% pq[cluster==cl | cluster==cl2 ]$contig] -> x
dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
x[, l := log10(0.1 + N)]
cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
morexfai[x, on="chr"] -> x
cen[x, on="chr"] -> x
x[, end_dist1 := pmin(pos1, len - pos1)]
x[, end_dist2 := pmin(pos2, len - pos2)]
x[, cen_dist1 := abs(cen_pos - pos1)]
x[, cen_dist2 := abs(cen_pos - pos2)]
x[, ldist := abs(pos1 - pos2)]
x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]
dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
y[, ctg1 := NULL]
prcomp(cor(y), scale=T, center=T)->pca
data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp1c1
pp1c1[, chr := j]
assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp1c1, on="contig"] -> pp1c1
pp1c1[, idx := 1:.N]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[cc<32, col := "red"]
pp1c1[cc>=32&cc<64,col:="blue"]
pp1c1[cc>=64&cc<96,col:="black"]
pp1c1[cc>=96,col:="orange"]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[, cluster := 5]
pdf(paste("haplotype_separation_HiC_pp",chrj,"c",cl,cl2,".pdf",sep=""), height=8, width=8)
plot(pp1c1$PC2, pp1c1$PC1, col=pp1c1$col, pch=20, las=1,bty='l', type="p", main=paste(chrj," cluster",cl,sep=""), ylab="PC1",xlab="PC2")
dev.off()
rbind(pq[!pp1c1$contig, on="contig"], pp1c1)-> pq1
print(sum(pq1[chr == chrj & cluster==1]$contig_length))
print(sum(pq1[chr == chrj & cluster==2]$contig_length))
print(sum(pq1[chr == chrj & cluster==3]$contig_length))
print(sum(pq1[chr == chrj & cluster==4]$contig_length))
print(sum(pq1[chr == chrj & cluster==5]$contig_length))
pp1c1[, idx := 1:.N]
pdf(paste("haplotype_separation_HiC_pp",chrj,"c",cl,cl2,"_repca.pdf",sep=""), height=4, width=10)
lapply(c("PC1", "PC2", "PC3", "PC4"), function(j){
 pp1c1[chr==chrj, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chrj), ylab=j,
            xlab="Hv syntenic position [Mb]")]
 pp1c1[chr==chrj, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(get(j), get(j))), by=idx]
})
dev.off()
return(pq1) 
}

Repcafor2check <- function(pq, chrj, cov,v, cl,cl2, prefix, plevel=0,pc=1,error_length=0)
{
pq=pq
j=chrj
dir.create(prefix)
lencl=sum(pq[chr == chrj & cluster==cl]$contig_length)
lencl2=sum(pq[chr == chrj & cluster==cl2]$contig_length)
setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
x[ctg1 %in%  pq[cluster==cl | cluster==cl2 ]$contig & ctg2 %in% pq[cluster==cl | cluster==cl2 ]$contig] -> x
dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
x[, l := log10(0.1 + N)]
cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
morexfai[x, on="chr"] -> x
cen[x, on="chr"] -> x
x[, end_dist1 := pmin(pos1, len - pos1)]
x[, end_dist2 := pmin(pos2, len - pos2)]
x[, cen_dist1 := abs(cen_pos - pos1)]
x[, cen_dist2 := abs(cen_pos - pos2)]
x[, ldist := abs(pos1 - pos2)]
x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]
dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
y[, ctg1 := NULL]
prcomp(cor(y), scale=T, center=T)->pca
data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp1c1
pp1c1[, chr := j]
assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp1c1, on="contig"] -> pp1c1
pp1c1[, idx := 1:.N]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[cc<32, col := "red"]
pp1c1[cc>=32&cc<64,col:="blue"]
pp1c1[cc>=64&cc<96,col:="black"]
pp1c1[cc>=96,col:="orange"]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
if(pc==1){
  pp1c1[PC1>plevel, cluster := cl]
  pp1c1[PC1<plevel, cluster := cl2]
}else if(pc==2){
  pp1c1[PC2>plevel, cluster := cl]
  pp1c1[PC2<plevel, cluster := cl2]
}else if(pc==3){
  pp1c1[PC3>plevel, cluster := cl]
  pp1c1[PC3<plevel, cluster := cl2]
}else if(pc==4){
  pp1c1[PC4>plevel, cluster := cl]
  pp1c1[PC4<plevel, cluster := cl2]
}
pdf(paste(prefix,"/haplotype_separation_HiC_pp",chrj,"c",cl,cl2,"_check.pdf",sep=""), height=8, width=8)
plot(pp1c1$PC2, pp1c1$PC1, col=pp1c1$col, pch=20, las=1,bty='l', type="p", main=paste(chrj," cluster",cl,sep=""), ylab="PC1",xlab="PC2")
dev.off()
pp1c1[, idx := 1:.N]
pp1c1->ppp
pq[, .(contig, clo=cluster)][ppp, on="contig"] -> ppp
pdf(paste(prefix,"/haplotype_separation_HiC_pp",chrj,"c",cl,cl2,"_repca_ckeck.pdf",sep=""), height=4, width=10)
lapply(c("PC1", "PC2", "PC3", "PC4"), function(j){
 ppp[chr==chrj, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chrj), ylab=j,
            xlab="Hv syntenic position [Mb]")]
 ppp[chr==chrj, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(get(j), get(j)),col=clo), by=idx]
 abline(h=plevel, col="blue")
})
dev.off()
rbind(pq[!pp1c1$contig, on="contig"], pp1c1)-> pq1
print(sum(pq1[chr == chrj & cluster==1]$contig_length))
print(sum(pq1[chr == chrj & cluster==2]$contig_length))
print(sum(pq1[chr == chrj & cluster==3]$contig_length))
print(sum(pq1[chr == chrj & cluster==4]$contig_length))
print(sum(pq1[chr == chrj & cluster==5]$contig_length))
checklencl=sum(pq1[chr == chrj & cluster==cl]$contig_length)
checklencl2=sum(pq1[chr == chrj & cluster==cl2]$contig_length)
if((abs(lencl - checklencl) <= error_length) | (abs(lencl - checklencl2) <= error_length)){
    print("Correct")
    print(abs(lencl - checklencl))
    print(abs(lencl - checklencl2))
}else{
    print("Need re-pca or check plevel")
}
}

####RE-PCA and check for 2 subgenome

readRDS(file="hparodii_assembly_v2_manual.Rds") -> pq1
cov[, .(contig=scaffold, cc)][pq1, on="contig"] -> pq1

##check chr1H phasing
Repcafor2check(pq1,"chr1H",cov,v,1,2,"chr1Hcheck",0)
Repcafor2check(pq1,"chr1H",cov,v,1,3,"chr1Hcheck",0)
Repcafor2check(pq1,"chr1H",cov,v,2,3,"chr1Hcheck",0)


####check chr2H phasing
Repcafor2check(pq1,"chr2H",cov,v,1,2,"chr2Hcheck",0)
Repcafor2check(pq1,"chr2H",cov,v,1,3,"chr2Hcheck",0)
Repcafor2check(pq1,"chr2H",cov,v,2,3,"chr2Hcheck",0)

####check chr3H phasing 
Repcafor2check(pq1,"chr3H",cov,v,1,2,"chr3Hcheck",0)
Repcafor2check(pq1,"chr3H",cov,v,1,3,"chr3Hcheck",0)
Repcafor2check(pq1,"chr3H",cov,v,2,3,"chr3Hcheck",0)

####check chr4H phasing 
Repcafor2check(pq1,"chr4H",cov,v,1,2,"chr4Hcheck",0)
Repcafor2check(pq1,"chr4H",cov,v,1,3,"chr4Hcheck",0.05)
Repcafor2check(pq1,"chr4H",cov,v,2,3,"chr4Hcheck",0.1)

####check chr5H phasing 
Repcafor2check(pq1,"chr5H",cov,v,1,2,"chr5Hcheck",0)
Repcafor2check(pq1,"chr5H",cov,v,1,3,"chr5Hcheck",0)
Repcafor2check(pq1,"chr5H",cov,v,2,3,"chr5Hcheck",0)

###check chr6H phasing 
Repcafor2check(pq1,"chr6H",cov,v,1,2,"chr6Hcheck",0)
Repcafor2check(pq1,"chr6H",cov,v,1,3,"chr6Hcheck",0)
Repcafor2check(pq1,"chr6H",cov,v,2,3,"chr6Hcheck",0.1)

####check chr7H phasing 
Repcafor2check(pq1,"chr7H",cov,v,1,2,"chr7Hcheck",0)
Repcafor2check(pq1,"chr7H",cov,v,1,3,"chr7Hcheck",0)
Repcafor2check(pq1,"chr7H",cov,v,2,3,"chr7Hcheck",0)

cal_chrlen(pq1,"chr1H",45)
cal_chrlen(pq1,"chr2H",45)
cal_chrlen(pq1,"chr3H",45)
cal_chrlen(pq1,"chr4H",45)
cal_chrlen(pq1,"chr5H",45)
cal_chrlen(pq1,"chr6H",45)
cal_chrlen(pq1,"chr7H",45)

saveRDS(pq1, file="hparodii_assembly_v2_haplotype_separation_v0.Rds")
