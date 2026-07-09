#!/bin/zsh

projectdir="/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hparodii"

reads=$projectdir/hifi_reads

prefix=${projectdir:r}/hparodii

find $reads | grep 'fastq.gz$' | xargs /hifiasm-0.15.5/hifiasm -o $prefix 2> ${prefix}_asm.err
