#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

#read chromosome lengths of Morex V3
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/200416_MorexV3_pseudomolecules.fasta.fai', sel=1:2, col.names=c("chr", "len"))->fai

#read centromere positions
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/MorexV3_centromere_positions.tsv')->cen
setnames(cen, c("chr", "cen_pos"))

##
# Import data and create assembly object
##

#read MorexV3 gene-based guide map ("pseudo-POPSEQ")
readRDS('/filer-dg/agruppen/seq_shared/mascher/hordeum_bulbosum_fb19_011_ccs_assembly_201111/Hb_FB19_011_hifiasm_201217/210129_MorexV3_HCgenes_pseudopopseq.Rds') -> pg

#read unitig lengths
f <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hjubatum_hifiasm.fasta.fai'
fread(f, head=F, select=1:2, col.names=c("scaffold", "length"))->fai

#read MorexV3 gene alignment (GMAP output) and merge with pseudo-POPSEQ table
fread('/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hic/hjubatum_hifiasm_MorexV3_HC_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))
fai[, .(scaffold, scaffold_length=length)][x[, .(css_contig=transcript, scaffold=contig, pos=start)], on="scaffold"] -> aln
pg[aln, on="css_contig"] -> aln

#read Hi-C map
dir <- '/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum/hic/hicreads'
fread(paste('find', dir, '| grep "fragment_pairs.tsv.gz$" | xargs zcat'),
      header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2"))->fpairs

#initialize assembly object, anchor to guide map and calculate physical coverage
init_assembly(fai=fai, cssaln=aln, fpairs=fpairs) -> assembly
anchor_scaffolds(assembly = assembly, popseq=pg, species="barley") -> assembly
add_hic_cov(assembly, binsize=1e4, binsize2=1e6, minNbin=50, innerDist=3e5, cores=20)->assembly

#save uncorrected assembly object
saveRDS(assembly, file="hjubatum_assembly.Rds")

##
# Check for and correct chimeras 
##

readRDS('hjubatum_assembly.Rds') -> assembly

#make diagnostic plot for one contig

#plot all > 3e6 contig 
assembly$info[length >= 3e6, .(scaffold, length)][order(-length)] -> ss

i=1
ss[i]$scaffold -> s 
assembly$cov[s, on='scaffold'] -> b

for (ni in 2:length(ss$scaffold)) {
i=ni
ss[i]$scaffold -> s 
rbind(b, assembly$cov[s, on='scaffold'])->b
} 

plot_chimeras(assembly=assembly, scaffolds=b,  species="barley", refname="Morex HC genes", autobreaks=F, mbscale=1, file="hjubatum_assembly_all_contig.pdf", cores=50)


##############################################################################################################################################################################################################


#set break point coordinates

i=29
ss[i]$scaffold -> s 
assembly$cov[s, on='scaffold'][bin >= 40e6 & bin <= 60e6][order(r)][1, .(scaffold, bin)] -> b

i=80
ss[i]$scaffold -> s 
rbind(b, assembly$cov[s, on='scaffold'][bin >= 15e6 & bin <= 20e6][order(r)][1, .(scaffold, bin)])->b




setnames(b, "bin", "br")
plot_chimeras(assembly=assembly, scaffolds=b, br=b, species="barley", refname="MorexV3 genes",  mbscale=1,
         file="hjubatum_assembly_chimeras_final.pdf", cores=30)

#implement the correction
break_scaffolds(b, assembly, prefix="contig_corrected_v1_", slop=1e4, cores=30, species="barley") -> assembly_v2

#save the object
saveRDS(assembly_v2, file="hjubatum_assembly_v2.Rds")

