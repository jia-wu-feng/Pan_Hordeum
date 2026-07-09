# Perform muscle alignment and trim it; input sequences include candidate gene and its best hits in references

module load muscle

for f in *all.fa ; do
	muscle -in $f -phyiout ${f/.fa/_muscle_in.phy}
	trimal -in ${f/.fa/_muscle_in.phy} -out trimmed_${f/.fa/_muscle_in.phy} -automated1
done
