# kir-genotyper
Algorithm for calling KIR gene alleles from BLASTed WXS or WGS reads.

# Installation
We recommend running KIR Genotyper using the Anaconda distribution of Python-3.4 or greater. Required python packages include numpy and pandas. For visualization purposes, matplotlib, seaborn, and sklearn are required.

```
git clone --recursive https://github.com/gaog94/KIRCLE
```

Then add the path to the kir-genotyper files to PYTHONPATH.

(((Include a part about setting up BLAST databases)))

# Description
KIR Genotyper is an algorithm that takes BLAST output from Whole Exome Sequencing (WXS) or Whole Genome Sequencing (WGS) reads mapping to a particular KIR gene and predicts the allelic genotype of that KIR gene in that sample.

As its input, KIR Genotyper uses nucleotide BLAST results in .csv format. Reads from WXS or WGS mapping to the KIR gene should be BLASTed against a database of all of that KIR gene's alleles.

The resulting read alignments are then fed into KIR Genotyper, which samples the reads with replacement and then runs an expectation maximization algorithm on each bootstrap until an arbitrary convergence threshold is attained. The resulting posterior probabilities are outputted and then averaged across bootstraps to determine final allele probability estimates.

If only a single allele exceeds a set probability threshold ``t``, then we call a homozygous solution.  If two  alleles exceed ``t``, then we call a heterozygous solution. If more than 2 alleles exceed ``t``, then we call the special null solution "No_Solution, No_Solution".

Output is written to file ``[input_filename]_genotype.txt``

# Usage

### Running Directly From Command Line
```
python process_single_cram.py fname.cram tag refLocations part nboot thresh alpha max_iter ncores
```
The following input parameters are use:
* fname.cram is the name of the input CRAM file.
* tag is the prefix of output files that KIRCLE will generate.
* refLocations is the path to the reference file containing KIR genomic coordinates in either hg19 or hg38.
* ``part`` is the number of proportion of reads bootstrapped with replacement with each iteration of the EM algorithm in KIRCLE.
* ``nboot`` is an integer representing the number of bootstraps to conduct. 
* ``thresh`` is the threshold used in the thresholding step.
* ``alpha`` is the value used to set the convergence criterion for the EM algorithm.
* ``max_iter`` is the maximum number of iterations the EM algorithm may run if convergence is not achieved.
* ``ncores`` is the number of cores KIRCLE should allocate.

Alternatively, if BLAST results are already available, KIRCLE may be run from the command line as follows:
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
