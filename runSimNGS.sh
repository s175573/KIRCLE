#!/bin/bash

kgene=KIR2DS1

python ../wessim/FastaCombos.py $kgene

# For each solution, simulate genomes and compute empirical probabilities
for f in $(ls */*.fa)
do
  echo "--> Working Through File " $f " <--"
  for i in `seq 100`;
  do
    simLibrary -n 2000 $f | simNGS s_1_4x.runfile.txt -p paired -o fasta >> ${f%.fa}_simreads_$i.fasta
  done
done
