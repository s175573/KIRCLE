#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 23:27:34 2019

@author: galengao
"""
import sys
import numpy as np
import pandas as pd

import EM_algorithm as em

def bootstrap(df, part=0.8):
    '''Given BLAST output as a dataframe, and a fraction of reads to bootstrap
    part, return a random bootstrapping of the BLASTed reads results.'''
    uniqueReads = df.qseqid.unique() # list of unique reads in BLAST results
    nSamp = int(part * len(uniqueReads)) # number of reads to sample
    if part < 1:
        bstrapReads = np.random.choice(uniqueReads, size=nSamp, replace=True)
    else:
        bstrapReads = np.random.choice(uniqueReads, size=nSamp, replace=True)

    return df[df.qseqid.isin(bstrapReads)]
    
def bootstrap_EM(df, part=0.8, n_boot=100, max_iter=100, alpha=0.000001):
    '''Given BLAST output as a dataframe, bootstrap fraction part of reads and
    run Expectation-Maximization algorithm on bootstrapped output n_boot times.
    Return mean of bootstrapped allele probability estimates & allele names.'''
    sortedAlleles = sorted(df.sseqid.unique())
    Pmat = np.zeros((n_boot, len(df.sseqid.unique())))
    for i in range(n_boot):
        df_bs = bootstrap(df, part=part)
        P, alleles = em.run_EM(df_bs, max_iter=max_iter, alpha=alpha)

        # careful: make sure alleles are in same order every time!!!
        sortedP = np.array([x for _,x in sorted(zip(alleles, P))])

        Pmat[i] = sortedP

    return np.mean(Pmat, axis=0), sortedAlleles


if __name__ == "__main__":
    fname = sys.argv[1]
#    fname = '../Workspace/KIR3DL2_blast.csv'
    
    df = pd.read_csv(fname, header=-1)
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',\
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = df[df.pident == 100]
    
    P, alleles = bootstrap_EM(df, part=0.8, n_boot=10, max_iter=100, alpha=0.000001)
    
    # write results to file
    outname = '_'.join(fname.split('_')[:-1]) + '_EMout.csv'
    pd.DataFrame(P, index=alleles, columns=['p']).to_csv(outname)
