module load bedtools

# Blast results (donors x Hordeum species and H. vulgare x Hordeum species) are compared to each other and filtered
# The filtration is based on donor locations ie. if a donor gene has a hit in both Hordeum (non-vulgare) and H. vulgare,
# then it's discarded (the overlap criterium is 20%)
# The remaining donor hits are cross-checked back with the blast results to get the subject (Hordeum species) coordinates

species_list=$1
donors=("Phal" "Pnotatum" "Pvag" "Sita")

for d in ${donors[@]} ; do

while read sp ; do

# from blast out6 output select donor locations - make a bed file (and remember to switch start-end if needed)
cut -f 1,7,8 $d*$sp*l150 | awk -v OFS='\t' '{if ($2>$3) print $1,$3,$2; else print $1,$2,$3}' | sort -k 1,1 -k 2,2n > ${d}_${sp}_hits_filtered.bed

# check with donor regions hit by vulgare - and leave only those NOT overlapping (overlap 20% required)
bedtools intersect -v -f 0.2 -a ${d}_${sp}_hits_filtered.bed -b ${d}_vulgare_hits_filtered.bed | sort -u > ${d}_${sp}_hits_without_Hv.bed

# see what hordeum regions that corresponds to
while read name start end ; do

awk -v n=$name -v s=$start -v e=$end '($1==n) && ((($7==s) && ($8==e)) || (($7==e) && ($8==s)))' $d*$sp*l150 > ${d}_${sp}_cds_hits_without_Hv

done < ${d}_${sp}_hits_without_Hv.bed

done < $species_list

done
