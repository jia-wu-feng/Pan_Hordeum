#!/bin/zsh

projectdir="/filer-dg/agruppen/dg6/fengj/panhordeum/assembly/hjubatum"

reads=$projectdir/hifi_reads

prefix=${projectdir:r}/hjubatum

find $reads | grep 'fastq.gz$' | xargs /filer-dg/agruppen/dg7/fengj/software/hifiasm-0.13/hifiasm -o $prefix 2> ${prefix}_asm.err
