# The identified HGT regions in Hordeum species (excl. vulgare) are cross-referenced with annotations that became available
# No overlap criterium (f parameter of bedtools intersect) is used, futher filtration will follow
# Results from all donors are pooled together, donor identification will be done downstream in phylogenetic analysis

module load samtools
module load bedtools

annot_dir=$1
species_list=$2
donors=("Phal" "Pnotatum" "Sita" "Pvag")

# do the intersections for each donor and each Hordeum species

while read species ; do

for d in ${donors[@]} ; do

# prepare a bed file
cut -f 2,9,10 ${d}_${species}_cds_hits_without_Hv | awk -v OFS='\t' '{if ($2>$3) print $1,$3,$2; else print $1,$2,$3}' | sort -k 1,1 -k 2,2n > ${d}_${species}_cds_hits_filtered.bed

# the actual job - intersect regions identified by blast with genes location
bedtools intersect -sorted -wo -a ${d}_${species}_cds_hits_filtered.bed -b ${annot_dir}/${species}*/${species}*.pgsb.r1.Mar2024.genes.bed > ${d}_${species}_cds_hits_filtered.annot_genes_no_f

# get only the names of genes (from annot)
cut -f 7 ${d}_${species}_cds_hits_filtered.annot_genes_no_f | sort -u > ${d}_${species}_cds_hits_filtered.annot_genes_no_f.names

done

# merge candidates from all donors
cat *${species}_cds_hits_filtered.annot_genes_no_f.names | sort -u > ${species}_candidates_no_f.names

# filter transposon-related and low confidence

while read line ; do

grep "$line" $annot_dir/${species}*/*csv | tr -d '"' | awk -F',' -v OFS='\t' '($3=="high") && ($4!="transposon-related") {print $1}'

done < ${species}_candidates_no_f.names > ${species}_candidates_no_f.transcript_names

done < ${species_list}
