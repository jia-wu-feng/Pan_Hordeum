module load newick-utils

# Perform a pre-screen of trees: if the neighbouring branch can be collapsed to one clade (and its PACMAD) then the case is clear,
# otherwise it needs to be checked manually

map=$1 # reference species - clade (PACMAD/BOP/outgroup) map for simplifying the tree
species_list=$2 # list of Hordeum species

while read species ; do

# after the renaming, tree will contain focus species in format H_pusillum
name=${species%_*}
short_name=${name/ordeum/}

for f in ${species}_trees/*tree ; do

# main op: reroot, then rename to basic groups, if neighbours have the same name, then collapse them to one branch, finally write out the name of sister clade
# ideally it will be one word, meaning that the whole sister clade was from the same group, any other names need to be reviewed
sister_clade=$(nw_reroot -l $f M_acuminat A_comosus_ | nw_rename - $map | nw_condense - | nw_clade -s - $short_name | sed 's/:.*//')

echo -e "$(basename $f)\t$sister_clade" >> ${species}_trees/overview.txt

done

done < $species_list

end
