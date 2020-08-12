#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 18:10:43 2019

@author: galengao
"""
import os
import sys
import time

import pandas as pd

from multiprocessing import Pool

import bootstrapper as bs
import genotyper as gt
import generate_multiple_KIR_bams as genKIR

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

def run_bootstrapper(scorefile, kgene, tag, part=0.8, nboot=100):
    '''Run bootstrapped EM algorithm on each BLAST score file to generate KIR
    allele probabilities. Return probability file.'''
    if os.stat(scorefile).st_size != 0: # make sure BLAST file is not empty
        P, alls = bs.bootstrap_BLAST_file(scorefile, kgene, pident=100, part=part, n_boot=nboot)
        df_out = pd.DataFrame(P, index=alls, columns=[tag])
        df_out.to_csv(scorefile[:-4]+'_calls.tsv', sep='\t')
        return df_out
    else:
        return pd.DataFrame([])


def run_genotyper(df, thresh=0.2):
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
    df_p = run_bootstrapper(scorefile, kgene, tag, part=part, nboot=nboot)
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



part = 0.25
nboot = 100
thresh = 0.2
ncores = 15
kgenes = ['KIR2DL1','KIR2DL2','KIR2DL3','KIR2DL4','KIR2DL5A','KIR2DL5B',\
          'KIR2DS1','KIR2DS2','KIR2DS3','KIR2DS4','KIR2DS5','KIR3DL1',\
          'KIR3DL2','KIR3DL3','KIR3DS1']
f = '../hg38_KIR_locations_noKV.tsv'

    
# get cram file and other runtime parameters
fname = sys.argv[1] # input CRAM file
tag = sys.argv[2] # output prefix
refLocations = sys.argv[3] # hg19 or hg38 KIR genomic coordinates
part = float(sys.argv[4]) # fraction of reads in each bootstrap
nboot = int(sys.argv[5]) # number of bootstraps
thresh = float(sys.argv[6]) # threshold parameter for genotyper
alpha = float(sys.argv[7]) # EM convergence criterion threshold
max_iter = int(sys.argv[8]) # maximum number of iterations of EM algorithm
ncores = int(sys.argv[9]) # number of cores to allocate

print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))

# split CRAM into multiple BAMs per KIR gene
start_time  = time.time()
print(' ||| Splitting CRAM file ||| ')
bamfiles = genKIR.split_master_bam(fname, tag, refLocations, hg='hg38')
print('CRAM file split')
print("--- %s seconds  ---"  % (time.time() - start_time))
print(' ')
print(' ')

# Serial processing alternative to parallel processing below
#outputs = [process_sample_multipleBAMs(b) for b in bamfiles]

# run town in parallel
print(' ||| Running Bootsrapper and Genotyper in Parallel... ||| ')
p = Pool(ncores)
outputs = p.map(process_sample_multipleBAMs, bamfiles)
dfs, sols = [x[0] for x in outputs], [x[1] for x in outputs]
print("--- %s seconds  ---"  % (time.time() - start_time))
print(' ')
print(' ')

# write KIR output to file
start_time  = time.time()
print(' ||| Writing outputs to file... ||| ')
write_KIR_outputs(tag, dfs, sols)
# write run parameters to file
write_param_outputs(tag, part, nboot, thresh)
print("--- %s seconds  ---"  % (time.time() - start_time))
print(' ')
print(' ')
