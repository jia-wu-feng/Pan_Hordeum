#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

#read assembly object. This is the hifiasm unitig assembly after breaking chimeras based on coverage drops
readRDS('hjubatum_assembly_v2.Rds') -> assembly_v2

##
# Positions contigs based on guide map
##

#read alignment of MorexV3 HC genes 
fread('./hjubatum_hifiasm_MorexV3_HC_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))

#read positions of MorexV3 genes
fread('/filer-dg/agruppen/seq_shared/mascher/morex_v3_annotation_200818/Hv_Morex.pgsb.Jul2020.HC_mrna_pos.txt')->p
setnames(p, c("chr", "chr_start", "chr_end", "transcript"))

#keep only genes aligned for >= 80 % of their length with 85 % identity, allow up to two alignment (one per haplotype), modify for tetraploids 
#####
##### important for pasing allotetraploids
#####
p[x, on="transcript"] -> px
px[alnlen >= 80 & id >= 85]->a
px[a[, .N, key=transcript][N <= 2]$transcript, on="transcript"]->aa

#convert original contig positions to corrected (assembly_v2) positions 
assembly_v2$info[, .(contig=orig_scaffold, start=orig_start, orig_start, scaffold)][aa, on=c("contig", "start"), roll=T] -> aa

#chromosome assignment
aa[, .N, key=.(chr, scaffold)][, p := N/sum(N), by=scaffold][order(-p)][!duplicated(scaffold)] -> cc

#keep contigs with at least 4 aligned genes, 75 % of aligned are from the major chromosome
cc[N >= 4 & p >= 0.75] -> cc

#check how much of the assembly can be assigned to chromosomes
assembly_v2$info[, .(scaffold, scaffold_length=length)][cc, on="scaffold"] -> cc

cc[scaffold_length>=3e5]->cc
sum(cc$scaffold_length)/1e6

####7383.251


#get the approximate chromosome positions (median of alignment coordinates)
aa[cc[, .(scaffold, chr)], on=c("chr", "scaffold")][, .(pos=median(as.numeric(chr_start)), pos_mad=mad(as.numeric(chr_start))), key=scaffold][cc, on="scaffold"] -> cc
cc[, mr := pos_mad / scaffold_length]  
setorder(cc, chr, pos) 
cc[, chr_idx := 1:.N, by=chr]
cc[, agp_pos := c(0, cumsum(scaffold_length[-.N])), by=chr]

#save results
saveRDS(cc, file="hjubatum_assembly_v2_Hv_guide.Rds")

##
# Positioning additional contigs by Hi-C
## 

#load guide map positions
readRDS('hjubatum_assembly_v2_Hv_guide.Rds') -> cc
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
saveRDS(a, file="hjubatum_assembly_v2_Hv_guide+HiClift.Rds")

#check how much of the assembly is placed
 sum(a[!is.na(chr)]$length)
#[1] 7466061148
 sum(a$length)
#[1] 7489204145
 7466061148/7489204145
#[1] 0.9911989

##
# Subgenome phasing
##

#read lengths of Morex chromosomes
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/200416_MorexV3_pseudomolecules.fasta.fai', sel=1:2, col.names=c("chr", "len"))->morexfai

#read centromere position in Morex 
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/MorexV3_centromere_positions.tsv')->cen
setnames(cen, c("chr", "cen_pos"))

#read scaffold posiitions and add them to the Hi-C link table
readRDS(file="hjubatum_assembly_v2_Hv_guide.Rds") -> cc
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
saveRDS(pp, file="hjubatum_assembly_v2_haplotype_separation_HiC.Rds")

#read coverage information and convert to correct assembly coordinates
fread(cmd="grep '^S' ../hjubatum.p_ctg.noseq.gfa | cut -f 2,5 | tr ':' '\\t' | cut -f 1,4", head=F) -> cov
setnames(cov, c("contig", "cc"))
cov[, contig := sub('l$', '', sub("ptg0*", "contig_", contig))]
assembly_v2$info[, .(contig=orig_scaffold, scaffold)][cov, on=c("contig")] -> cov
setorder(cc[, .(scaffold, scaffold_length, chr, pos)][cov, on="scaffold"], chr, pos) -> cov
cov[, idx := 1:.N]
saveRDS(file="hjubatum_assembly_v2_cov.Rds", cov)

pp-> pq

#manually define cuts between haplotype for each chromosome
#2x contigs get "haplotype 3", i.e. present in both 1 and 2
cov[, .(contig=scaffold, cc)][pq, on="contig"] -> pq
pq[chr == "chr1H" & PC1 > 0, hap := 1]
pq[chr == "chr1H" & PC1 <= 0, hap := 2]
pq[chr == "chr2H" & PC2 > 0.04, hap := 1]
pq[chr == "chr2H" & PC2<= 0.04, hap := 2]
pq[chr == "chr3H" & PC1 > 0, hap := 1]
pq[chr == "chr3H" & PC1 <= 0, hap := 2]
pq[chr == "chr4H" & PC1 > 0, hap := 1]
pq[chr == "chr4H" & PC1 <= 0, hap := 2]
pq[chr == "chr5H" & PC1 > 0, hap := 1]
pq[chr == "chr5H" & PC1 <= 0, hap := 2]
pq[chr == "chr6H" & PC1 > 0, hap := 1]
pq[chr == "chr6H" & PC1 <= 0, hap := 2]
pq[chr == "chr7H" & PC1 > 0, hap := 1]
pq[chr == "chr7H" & PC1 <= 0, hap := 2]
pq[, col := "black"]
pq[hap == 1, col := "red"]
pq[hap == 2, col := "blue"]



fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/200416_MorexV3_pseudomolecules.fasta.fai', sel=1:2, col.names=c("chr", "len"))->morexfai

#plot results (PC1 score and coverage)
pdf("210210_assembly_v2_haplotype_separation_HiC.pdf", height=8, width=10)
par(mfrow=c(2, 1))
lapply(c("PC1"), function(j){
 lapply(morexfai[1:7, chr], function(i){
  pq[chr == i, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chr[1]), ylab="PC1",
        xlab="Hv syntenic position [Mb]", xlim=c(0, morexfai[i, len/1e6, on="chr"]))]
   pq[chr == i, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(PC1, PC1), col=col), by=idx]
  abline(v=c(0, morexfai[i, len/1e6, on="chr"]), col="blue")
  cov[chr == i, plot(pos/1e6, las=1, bty='l', cc, type='n', main=sub("chr", "", chr[1]), ylab="coverage",
        xlab="Hv syntenic position [Mb]", xlim=c(0, morexfai[i, len/1e6, on="chr"]))]
   cov[chr == i, lines(lwd=3, c(pos/1e6, (pos + scaffold_length)/1e6), c(cc, cc), col=1), by=idx]
 })
})

dev.off()


pdf("hjubatum_separation_HiC_all.pdf", height=4, width=10)
par(mfrow=c(2,4))
lapply(c("PC1", "PC2", "PC3", "PC4"), function(j){
 lapply(morexfai[1:7, chr], function(i){
  pq[chr == i, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chr[1]), ylab=j,
        xlab="Hv syntenic position [Mb]")]
  pq[chr == i, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(get(j), get(j))), by=idx]
 })
 plot(axes=F, xlab="", ylab="", 0, type='n')
})
dev.off()

pdf("hjubatum_separation_HiC_all_pc1_2.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pq[chr == i, plot(PC2, PC1, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
 })

dev.off()


saveRDS(pq, file="hjubatum_assembly_v2_haplotype_separation_v0.Rds")




