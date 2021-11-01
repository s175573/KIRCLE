#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 23:27:34 2019

@author: galengao
"""
from __future__ import division

import numpy as np
import pandas as pd
import EM_algorithm as em

from collections import Counter


def bootstrap(df, part=0.5):
    '''Given BLAST output as a dataframe, and a fraction of reads to bootstrap
    part, return a random bootstrapping of the BLASTed reads results.'''
    uniqueReads = df.qseqid.unique() # list of unique reads in BLAST results
    nSamp = int(part * len(uniqueReads)) + 1 # #reads in sample -- always >= 1
    bstrapReads = np.random.choice(uniqueReads, size=nSamp, replace=True)

    bsreadcount = Counter(bstrapReads)

    # handle duplicate reads
    dfs = []
    for i in range(max(bsreadcount.values())):
        r = [x for x in bsreadcount if bsreadcount[x] >= i+1]
        df_r = df[df.qseqid.isin(r)]
        df_r.qseqid = df_r.qseqid + '_' + str(i)
        dfs.append(df_r)

    return pd.concat(dfs)

def bootstrap_EM(df, kgene, part=0.5, n_boot=100, maxIter=1000, alpha=1e-5):
    '''Given BLAST output as a dataframe, bootstrap fraction part of reads and
    run Expectation-Maximization algorithm on bootstrapped output n_boot times.
    Return mean of bootstrapped allele probability estimates & allele names.'''
    # read in list of possible alleles from file, not from BLAST outtput
    df_tr = pd.read_csv('../ref_files/KIRallele_translation.txt', sep='\t', index_col=0)
    tdict = df_tr['Colon_Name'].to_dict()

    # load list of alleles and possible 2-allele combinations
    df_kirvars = pd.read_csv('../ref_files/KIR_variants.txt', sep='\t')
    sortedAlleles = sorted([tdict[a] for a in df_kirvars[kgene].dropna()])

    # set up bootstrapping
    Pmat = np.zeros((n_boot, len(sortedAlleles)))
    for i in range(n_boot):
        df_bs = bootstrap(df, part=part)

        if len(df_bs) != 0:
            P, alleles = em.run_EM(df_bs, maxIter=maxIter, alpha=alpha)

            pdict = {a:P[i] for i, a in enumerate(alleles)}
            for a in sortedAlleles:
                if a not in pdict:
                    pdict[a] = 0

            # Make sure alleles are in same order every time
            sortedP = np.array([pdict[a] for a in sortedAlleles])

            if len(sortedP) == len(Pmat[i]):
                Pmat[i] = sortedP

    return np.mean(Pmat, axis=0), sortedAlleles

def bootstrap_BLAST_file(fname, kgene, pident=100, part=0.5, n_boot=100, maxIter=1000, alpha=1e-5):
    '''Given file of BLAST scores (output format 10), run bootstrap_EM() to
    compute mean of bootstrapped allele probability estimates & allele names.'''
    df = pd.read_csv(fname, header=None)
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',\
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = df[df.pident == pident]

    # double check to make sure BLAST results aren't empty
    if len(df != 0):
        P, alls = bootstrap_EM(df, kgene, part=part, n_boot=n_boot, maxIter=maxIter, alpha=alpha)
    else: # otherwise generate null results
        P, alls = np.empty(0), []

    return P, alls