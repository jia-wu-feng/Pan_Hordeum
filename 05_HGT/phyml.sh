#!/bin/bash

# Construct a phylogenetic tree from the alignment

module add phyml-3.1

for f in *_in.phy ; do

phyml -i $f -m GTR -b 1000

done
