#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 06:23:16 2019

@author: galengao
"""
import sys
import numpy as np
import pandas as pd

import multiprocessing as mp

import bootstrapper as bs
import genotyper

def process_multiple_thresholds(f):
    print('Characterizing file ' + f + ' ...')
    df = pd.read_csv(f, header=-1)
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',\
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = df[df.pident == 100]

    P, alleles = bs.bootstrap_EM(df, part=0.05, n_boot=100, max_iter=100, alpha=0.000001)
    pd.DataFrame(P, index=alleles).to_csv('_'.join(f.split('_')[:-1]) + '_belief.csv')

    ts = np.linspace(0.02, 0.4, 20)
    tResults = []
    for t in ts:
        a1, a2 = genotyper.genotype_bootstrap(P, alleles, thresh=t)
        tResults.extend([a1, a2])

    return tResults    

fnames = sys.argv[1:]


pool = mp.Pool(processes=32)
results = pool.map(process_multiple_thresholds, fnames)
    

print('Writing Results to File...')
ts = np.linspace(0.02, 0.4, 20)
cols = np.array([["{0:.2f}".format(t)+'_a1', "{0:.2f}".format(t)+'_a2'] \
                  for t in ts]).flatten()
df_out = pd.DataFrame(results, columns=cols)
print(df_out.head())
directory = '/'.join(fnames[0].split('/')[:-1])
df_out.to_csv(directory+'Genotyper_results_different_thresholds.csv')