#!/bin/bash

### This script uses simNGS to generate 1000 simulated genomes and then runs
### an expectation maximization algorithm to genotype each simulated genome

run_simulations () {
  sol=$(dirname $1)
  n=$2
  nproc=$3

  # Generate n simulated genomes using simNGS; parallelized
  echo " ||| Generating " $n " Genome Simulations... ||| "
  seq $n | parallel -j $nproc "simLibrary -n 4000 $f | simNGS s_1_4x.runfile.txt \
    -p paired -o fasta >> ${f%.fa}_simreads_{}.fasta"

  # BLAST n simulated genomes
  db=${kgene}_db
  echo " ||| BLASTing " $n " Simulated Genomes... ||| "
  ls $sol/*.fasta | parallel -j $nproc "blastn -db $db -outfmt 10 -max_hsps 50 -query {} > {.}_blast.csv"

  rm $sol/*.fasta

  # Run Genotyper algorithm to estimate our solution
  echo " ||| Running Genotyper on " $n " BLAST results... ||| "
  #ls $sol/*.csv | parallel -j $nproc python genotyper.py {} 100
  python genotyper_simulations.py $sol/*_blast.csv
}

kgene=$1

# Build directory structure with one genotype solution per directory
python FastaCombos.py $kgene

# For each solution, simulate genomes and compute empirical probabilities
for f in $(ls $kgene*/*.fa)
do
  echo "--> Working Through File " $f " <--"
  run_simulations $f 500 32
done
