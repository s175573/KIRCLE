#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 18:10:43 2019

@author: galengao
"""
import os
import time
import argparse

import pandas as pd

from multiprocessing import Pool

import bootstrapper as bs
import genotyper as gt
import make_KIR_minibams as genKIR

def convert_to_fasta(bamfile):
    '''Unpack reads from BAM file (xxx.bam) into FASTA file (xxx.fasta).
    Returns name of FASTA file.'''
    stem = bamfile[:-4]
    os.system('samtools fasta  ' + bamfile + ' > ' + stem + '.fasta')

    return stem+'.fasta'

def run_blast(fastafile, kgene):
    '''BLAST reads from inputted BAMfile against a database of all reference
    alleles for the inputted KIR gene; outputted as a .csv. Returns CSV name.'''
    kDB = kgene + '_db'
    outfile = fastafile[:-6]  + '.csv'
    cmd = 'blastn -db ' +kDB+ ' -outfmt 10 -query ' +fastafile+ ' > ' +outfile
    os.system(cmd)

    return outfile

def run_bootstrapper(scorefile, kgene, tag, part=0.5, nboot=100, alpha=1e-05, maxIter=1000):
    '''Run bootstrapped EM algorithm on each BLAST score file to generate KIR
    allele probabilities. Return probability file.'''
    if os.stat(scorefile).st_size != 0: # make sure BLAST file is not empty
        P, alls = bs.bootstrap_BLAST_file(scorefile, kgene, pident=100, part=part, \
                                          n_boot=nboot, alpha=alpha, maxIter=maxIter)
        df_out = pd.DataFrame(P, index=alls, columns=[tag])
        df_out.to_csv(scorefile[:-4]+'_calls.tsv', sep='\t')
        return df_out
    else:
        return pd.DataFrame([])


def run_genotyper(df, thresh=0.25):
    '''Run KIR genotyper on each BLAST score file.'''
    if len(df) != 0: # make sure allele probability file is not empty
        sol = gt.genotype_bootstrap(df.values.ravel(), df.index, thresh=thresh)
    else:
        sol = ('No Solution', 'No Solution')

    return sol

def process_sample_multipleBAMs(bamfile):
    '''Given inputted BAM file, run bootstrapper & genotyper. Individual
    functions write output to disc.'''
    # extract the KIR gene we're dealing with from the bamfile name
    stem = bamfile[:-4]
    kgene = stem.split('_')[-1]
    tag = stem.split('_')[0]

    print(' >>> ' + kgene + ' <<< ')
    # run pipeline on the BAMfile
    fastafile = convert_to_fasta(bamfile)
    scorefile = run_blast(fastafile, kgene)
    df_p = run_bootstrapper(scorefile, kgene, tag, part=part, nboot=nboot, \
                            alpha=alpha, maxIter=maxIter)
    sol = run_genotyper(df_p, thresh=thresh)

    return df_p, sol

def write_KIR_outputs(tag, df_outs, sols):
    '''Write outputs to disc: [tag]_probabilities.txt and [tag]_genotypes.txt'''
    # write aggregated allele probabilities to file
    if len(df_outs) != 0:
        df_out = pd.concat(df_outs)
    else:
        df_out = pd.DataFrame([])
    df_out = df_out.loc[~df_out.index.duplicated(keep='first')]
    df_out.to_csv(tag+'_probabilities.txt', sep='\t')

    # write aggregated genotype solutions to file
    outname = tag+'_genotypes.txt'
    with open(outname, 'w') as tsv:
        for k, s in zip(kgenes, sols):
            tsv.write(k+ '\t' + s[0] + '\t' + s[1] + '\n')

def write_param_outputs(tag, part, nboot, thresh):
    '''Write runtime parameters to disc: [tag]_runParams.txt'''
    outname = tag+'_runParams.txt'
    with open(outname, 'w') as tsv:
        for n,p in zip(['part', 'nboot', 'thresh'], [part, nboot, thresh]):
            tsv.write(n + '\t' + str(p) + '\n')


kgenes = ['KIR2DL1','KIR2DL2','KIR2DL3','KIR2DL4','KIR2DL5A','KIR2DL5B',\
          'KIR2DS1','KIR2DS2','KIR2DS3','KIR2DS4','KIR2DS5','KIR3DL1',\
          'KIR3DL2','KIR3DL3','KIR3DS1']


# parse runtime parameters, store as variables, and echo parameter values
purpose = "Perform KIR alelle inference on a single CRAM or BAM using KIRCLE"
parser = argparse.ArgumentParser(description=purpose)

parser.add_argument('-i', '--input', type=str, help="CRAM or BAM input file")
parser.add_argument('-o', '--tag', type=str, \
                    help="Prefix to append to all output files")
parser.add_argument('-l', '--loci', type=str, \
                    default="../ref_files/hg38_KIR_locations_noKV.tsv", \
                    help="File of KIR genomic coordinates/loci")
parser.add_argument('-g', '--genome', choices=['hg19','hg38'], default='hg38', \
                    help="Reference genome that input BAM/CRAM is aligned to")
parser.add_argument('-p', '--partition', type=float, default=0.5, \
                    help="Fraction of reads to include in each bootstrap")
parser.add_argument('-b', '--bootstraps', type=int, default=100, \
                    help="Number of bootstraps to perform")
parser.add_argument('-t', '--threshold', type=float, default=0.25, \
                    help="Threshold parameter for Genotyper algorithm")
parser.add_argument('-a', '--alpha', type=float, default=0.00001, \
                    help="Threshold for EM algorithm convergence")
parser.add_argument('-m', '--maxIterations', type=int, default=1000, \
                    help="Maximum # of iterations to run EM algorithm")
parser.add_argument('-c', '--cores', type=int, default=1, \
                    help="Number of cores to allocate for parallel processing")
args = parser.parse_args()

fname = args.input
tag = args.tag
refLocations = args.loci
hg = args.genome
part = args.partition
nboot = args.bootstraps
thresh = args.threshold
alpha = args.alpha
maxIter = args.maxIterations
ncores = args.cores

print("Input BAM/CRAM file: " + fname)
print("Output prefix tag: " + tag)
print("Reference file of KIR loci: " + refLocations)
print("Human genome reference: " + hg)
print(f"Fraction of reads in each bootstrap: {part}")
print(f"Number of bootstraps to perform: {nboot}")
print(f"Genotyper threshold parameter: {thresh}")
print(f"Threshold for EM algorithm convergence: {alpha}")
print(f"Maximum number of iterations to run EM algorithm: {maxIter}")
print(f"Number of cores to allocate for parallel processing: {ncores}")
print(' ')
print(' ')


# split CRAM into multiple BAMs per KIR gene
start_time  = time.time()
print(' ||| Splitting input file ||| ')
bamfiles = genKIR.split_master_bam(fname, tag, refLocations, hg=hg)
print('Input file split')
print("--- %s seconds  ---"  % (time.time() - start_time))
print(' ')
print(' ')

# Serial processing alternative to parallel processing below
#outputs = [process_sample_multipleBAMs(b) for b in bamfiles]

# Run KIRCLE in parallel
print(' ||| Running Bootsrapper and Genotyper in Parallel... ||| ')
p = Pool(ncores)
outputs = p.map(process_sample_multipleBAMs, bamfiles)
dfs, sols = [x[0] for x in outputs], [x[1] for x in outputs]
print("--- %s seconds  ---"  % (time.time() - start_time))
print(' ')
print(' ')

# write KIR outputs to file
start_time  = time.time()
print(' ||| Writing outputs to file... ||| ')
write_KIR_outputs(tag, dfs, sols)
# write run parameters to file
write_param_outputs(tag, part, nboot, thresh)
print("--- %s seconds  ---"  % (time.time() - start_time))
print(' ')
print(' ')