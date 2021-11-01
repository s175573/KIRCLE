#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  4 00:03:47 2019

@author: galengao
"""
import numpy as np

import bootstrapper as bs

def genotype_bootstrap(P, alleles, thresh=0.25):
    '''Given bootstrap results, use theshold thresh to determine specific
    solution (homo, hetero, indeterminate).'''
    # if we only have a single allele surpassing the threshold: homozygous
    if sum(P > thresh) == 1:
        mlv = alleles[np.argmax(P)]
        return (mlv, mlv)
    # otherwise, call a heterozygous solution, if only 2 between thresh & 1-thresh
    elif sum(P > thresh) == 2:
        hetAlleles = np.where(P > thresh, alleles, '')
        v1, v2 = [x for x in hetAlleles if x != '']
        return (v1, v2)
    # if >2 variants at het level or nothing at that level, then we're confused
    else:
        return ('No Solution', 'No Solution')

def genotype_results(df, kgene, thresh=0.25, part=0.5, n_boot=100, max_iter=1000):
    '''Given BLAST results, use threshold thresh to determine if we have a
    homozygous, heterozygous, or indeterminate solution.'''
    # Ensure columns named correctly; only consider 100% matched reads
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',\
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = df[df.pident == 100]

    # Run bootstrap EM algorithm on BLASTresults
    P, alleles = bs.bootstrap_EM(df, kgene, part=part, n_boot=n_boot, \
                                 max_iter=max_iter, alpha=0.000001)
    return genotype_bootstrap(P, alleles, thresh=thresh)