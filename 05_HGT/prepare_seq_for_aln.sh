source ~/toolbox.sh
start

module load bedtools

### requires references already indexed by samtools

while read species ; do

if [ ${species:0:1} == "#" ] ; then continue ; fi

echo $species

base=/storage/plzen1/home/chudma04/hordeum_panicoid_hgt/phylo_control/$species
blast_out=$base/blast_refs_stage_II
workdir=$base/workdir_extract_seq_stage_II
refdir=/storage/plzen1/home/chudma04/hordeum_panicoid_hgt/phylo_control/refs_stage_II
muscle=$base/muscle_aln_stage_II
candidates=$base/representative_candidates.fa

mkdir -p $workdir $muscle

while read name ; do

for out in $blast_out/candidates*out6 ; do

	part_ref=${out/*candidates_vs_/}
	ref=${part_ref/.out6/}

	ref_loc=$(grep "$name" $out | head -n 1 | cut -f 2)
	if [ ! -z "$ref_loc" ] ; then
		grep $ref_loc $refdir/$ref.fa.fai | awk -v r="$ref" -v OFS='\t' '{print $1,0,$2,r}' > $workdir/${name}_${ref}.bed
		bedtools getfasta -fi $refdir/$ref.fa -bed $workdir/${name}_${ref}.bed -fo $workdir/${name}_${ref}.fa -name
	fi

done

grep -A 1 "$name" $candidates | sed 's/'$name'/'${species%_*}'/' > $workdir/${name}.fa
cat $workdir/${name}*.fa > $muscle/${name}_all.fa

done < $base/only_candidates_in_orthogroups

done < species_list_missing

end
