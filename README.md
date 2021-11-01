# KIRCLE
Algorithm for calling KIR alleles inference from Whole Exome Sequencing data.

# Installation
We recommend running KIRCLE using Python-3.4 or greater. Required python packages include collections, numpy, pandas, pysam, multiprocessing, and warnings. To install KIRCLE, run the following code:

```
git clone --recursive https://github.com/gaog94/KIRCLE
```

Then add the pulled directory to PYTHONPATH.

To construct KIR-gene-specific BLAST databases, first download KIR nucleotide sequences in FASTA format from EBI (https://www.ebi.ac.uk/ipd/kir/download.html). Separate the downloaded file into separate fasta files each containing variants corresponding to a single gene for each KIR gene.

Then construct a BLAST database for each KIR gene by running:
```
Makeblastdb -in KIRgene.fasta.txt -parse_seqids -dbtype nucl -out KIRgene_db
```

where KIRgene.fasta.txt represents the fasta file containing the KIR-gene-specific sequences to be contained in the database, and KIRgene_db represents the name of the resulting BLAST database (e.g. "KIR2DL4_db").

BLAST may be obtained at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/.

# Description
KIRCLE is an algorithm that takes Whole Exome Sequencing (WXS) or Whole Genome Sequencing (WGS) data and infers the sample's KIR allele genotypes.

KIRCLE takes either an aligned SAM, BAM, or CRAM file as its input. It then performs nucleotide BLAST against a database of known KIR allele variants to generate BLAST alignment scores.

The resulting BLAST alignment scores are then fed into the bootstrapper algorithm, which samples the reads with replacement and then runs an expectation maximization algorithm on each bootstrap until a heuristically selected convergence threshold (``alpha``) is attained. The resulting posterior probabilities are outputted and then averaged across bootstraps to determine final allele probability estimates (outputted as ``[tag]_probabilities.txt``).

If only a single allele exceeds a set probability threshold ``threshold``, then we call a homozygous solution.  If two  alleles exceed ``threshold``, then we call a heterozygous solution. If more than 2 alleles exceed ``threshold``, then we call the special null solution "No_Solution, No_Solution".

Output is written to file ``[tag]_genotypes.txt``

<!---
# Diagnostics
A diagnostics script has been packaged with this algorithm. This script ensures proper functioning of the EM algorithm, bootstrapper, and the thresholding algorithm. To run diagnostics on KIR Genotyper:
```
python diagnostics_package.py
```
--->

# Usage

### Running KIRCLE Directly From Command Line
```
python KIRCLE.py -i filename -o tag -l refLoci -g hgVersion -p partition -b bootstrap -t threshold -a alpha -m maxIter -c nCores
```
The following input parameters are used:
* ``-i`` ``filename`` is the name of the input SAM/BAM/CRAM file.
* ``-o`` ``tag`` is the prefix of output files that KIRCLE will generate.
* ``-l`` ``refLoci`` is the path to the reference file containing KIR genomic coordinates/loci. (default "../ref_files/hg38_KIR_locations_noKV.tsv")
* ``-g`` ``hgVersion`` is either "hg19" or "hg38" depending on which reference the input file is aligned against. (default "hg38")
* ``-p`` ``partition`` is the proportion of reads bootstrapped with replacement in each iteration of the EM algorithm in KIRCLE. (default 0.5)
* ``-b`` ``bootstrap`` is an integer representing the number of bootstraps to conduct. (default 100)
* ``-t`` ``threshold`` is the threshold used in the thresholding step. (default 0.25)
* ``-a`` ``alpha`` is the threshold to use to determine convergence of the EM algorithm. (default 1e-05)
* ``-m`` ``maxIter`` is the maximum number of iterations the EM algorithm may run before reporting that convergence was not achieved. (default 1000)
* ``-c`` ``nCores`` is the number of cores KIRCLE should allocate for parallel processing. (default 1)

# Attributions
KIRCLE is developed and maintained by Galen Gao (galen [period] gao [at] utsw [period] edu) -- Bo Li lab -- UT Southwestern, Dallas, TX.
