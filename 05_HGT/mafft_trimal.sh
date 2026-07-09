# Perform alignment and its trimming using candidate sequences and best hits from each reference as input
# Pre-filtering of alignments: hits from at least 5 references, including 1 of panicoid species, needed

module load mafft

mafft=$1 # dir with input fasta sequences, will be also used as output dir
pan=$2 # list of panicoid species names

cd $mafft

for f in *all.fa ; do
        # references needed
        ref_count=$(grep -c ">[^H]" $f)
        if [ $ref_count -lt 5 ] ; then continue ; fi


	# check for panicoid reference
	panicoid_count=$(grep -c -f $pan $f)
	if [ $panicoid_count -lt 1 ] ; then continue ; fi

	mafft --auto --phylipout $f > ${f/.fa/_mafft_in.phy}
        trimal -in ${f/.fa/_mafft_in.phy} -out trimmed_${f/.fa/_mafft_in.phy} -automated1

done 

cd -
