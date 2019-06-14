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

def bootstrap_readnames(uniqueReads, part=0.8):
    '''Given a list of read names, pick a few to return.'''
    nSamp = int(part * len(uniqueReads)) + 1 # n_reads to sample -- always >= 1
    return np.random.choice(uniqueReads, size=nSamp, replace=True)

def bootstrap(df, part=0.8):
    '''Given BLAST output as a dataframe, and a fraction of reads to bootstrap
    part, return a random bootstrapping of the BLASTed reads results.'''
    uniqueReads = df.qseqid.unique() # list of unique reads in BLAST results
    nSamp = int(part * len(uniqueReads)) + 1 # # reads to sample -- always >= 1
    if part < 1:
        bstrapReads = np.random.choice(uniqueReads, size=nSamp, replace=True)
    else:
        bstrapReads = np.random.choice(uniqueReads, size=nSamp, replace=True)

    return df[df.qseqid.isin(bstrapReads)]
    
def bootstrap_EM(df, kgene, part=0.8, n_boot=100, max_iter=100, alpha=0.000001):
    '''Given BLAST output as a dataframe, bootstrap fraction part of reads and
    run Expectation-Maximization algorithm on bootstrapped output n_boot times.
    Return mean of bootstrapped allele probability estimates & allele names.'''
#    sortedAlleles = sorted(df.sseqid.unique())
    # read in list of possible alleles from file, not from BLAST outtput
    df_tr = pd.read_csv('../KIRallele_translation.txt', sep='\t', index_col=0)
    tdict = df_tr['Colon_Name'].to_dict()
    
    # load list of alleles and 2-allele combos
    df_kirvars = pd.read_csv('../KIR_variants.txt', sep='\t')
    sortedAlleles = sorted([tdict[a] for a in df_kirvars[kgene].dropna()])

    # set up bootstrapping    
    uniqueReads = df.qseqid.unique()
    alignDict = {r:df[df.qseqid == r] for r in uniqueReads}        
    
    Pmat = np.zeros((n_boot, len(sortedAlleles)))
    for i in range(n_boot):
        bsreads = bootstrap_readnames(uniqueReads)
        df_bs = pd.concat([alignDict[r] for r in bsreads], ignore_index=True) # make sure it's w/replace
#        df_bs = bootstrap(df, part=part)

        if len(df_bs) != 0:
            P, alleles = em.run_EM(df_bs, max_iter=max_iter, alpha=alpha)
            
            pdict = {a:P[i] for i, a in enumerate(alleles)}
            for a in sortedAlleles:
                if a not in pdict:
                    pdict[a] = 0
    
            # careful: make sure alleles are in same order every time!!!
#            sortedP = np.array([x for _,x in sorted(zip(alleles, P))])
            sortedP = np.array([pdict[a] for a in sortedAlleles])
    
            if len(sortedP) == len(Pmat[i]):
                Pmat[i] = sortedP

    return np.mean(Pmat, axis=0), sortedAlleles


if __name__ == "__main__":
#    fname = sys.argv[1]
    fname = '../Workspace/KIR3DL1_blast.csv'
    fname = 'KIR2DS2_GenotyperResults/KIR2DS2*0010101_KIR2DS2*0010108_simreads_42_blast.csv'

    
    df = pd.read_csv(fname, header=-1)
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',\
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = df[df.pident == 100]
    
    P, alleles = bootstrap_EM(df, part=0.3, n_boot=100, max_iter=100, alpha=0.000001)
    
    # write results to file
#    outname = '_'.join(fname.split('_')[:-1]) + '_EMout.csv'
#    pd.DataFrame(P, index=alleles, columns=['p']).to_csv(outname)
    
    import matplotlib.pyplot as plt
    # plot EM algorithm output in graphically appealing form
    plt.plot(P, 'bo')
    plt.xticks(range(len(alleles)), alleles, rotation=30, ha='right')
    plt.xlabel('Variant')
    plt.ylabel('Posterior Probability')
#    plt.title('TCGA-OR-A5J2-10A; KIR3DL2', fontsize=20)
    plt.title('KIR2DS2_Heterozygous_Simulation', fontsize=20)
