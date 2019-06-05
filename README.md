# kir-genotyper
Algorithm for calling KIR gene alleles from BLASTed WXS or WGS reads.

# Installation
We recommend running KIR Genotyper using the Anaconda distribution of Python-3.4 or greater. Required python packages include numpy and pandas. For visualization purposes, matplotlib, seaborn, and sklearn are required.

```
git clone --recursive https://github.com/gaog94/kir-genotyper
```

Then add the path to the kir-genotyper files to the path.

# Description
KIR Genotyper is an algorithm that takes BLAST output from Whole Exome Sequencing (WXS) or Whole Genome Sequencing (WGS) reads mapping to a particular KIR gene and predicts the allelic genotype of that KIR gene in that sample.

As its input, KIR Genotyper uses nucleotide BLAST results in .csv format. Reads from WXS or WGS mapping to the KIR gene should be BLASTed against a database of all of that KIR gene's alleles.

The resulting read alignments are then fed into KIR Genotyper, which samples the reads with replacement and then runs an expectation maximization algorithm on each bootstrap until an arbitrary convergence threshold is attained. The resulting posterior probabilities are outputted and then averaged across bootstraps to determine final allele probability estimates.

If only a single allele exceeds a set probability threshold ``t``, then we call a homozygous solution.  If two  alleles exceed ``t``, then we call a heterozygous solution. If more than 2 alleles exceed ``t``, then we call the special null solution "No_Solution, No_Solution".

Output is written to file ``[input_filename]_genotype.txt``

# Usage

### Running Directly From Command Line
```
python genotyper.py BLAST_output.csv <n_boot>
```
BLAST_output.csv is the result of running nucleotide blast on WXS or WGS reads mapping to a particular KIR gene against a database of alleles corresponding to that KIR gene. ``<n_boot>`` is an integer representing the number of bootstraps to conduct.

### Running as a Python Package
```
import genotyper as gt

# Load BLAST results in .csv format
df = pd.read_csv(inputFileName, header=-1)

# Run KIR genotyper
a1, a2 = gt.genotype_results(df, thresh=0.02, part=0.8, n_boot=20, max_iter=100)
```
``genotyper`` can be imported as a python package so that its constituent functions may be called within a python script.

# Attributions
This kir genotyper is developed and maintained by Galen Gao (galen.gao@utsw.edu) -- Bo Li lab -- UT Southwestern, Dallas, TX.
