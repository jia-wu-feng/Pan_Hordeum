blast_query=input_species_cds # potential donor (P.hallii, Pa. notatum, S. italica
ref_file=hordeum_list # list of Hordeum species
 

mapfile -t list1 < "$blast_query"  # Read lines from input_species_cds
mapfile -t list2 < "$ref_file"    # Read lines from hordeum_list

for f1 in "${list1[@]}"
do
    for f2 in "${list2[@]}"
    do
        echo "Running BLAST for: $f1 and ${f2}.high.cds.fa"
        blastn -query "$f1" -db "${f2}.cds.fa" -outfmt 6 -out blast_out/$output_file -num_threads 10 -evalue 1e-10 
    done
done
