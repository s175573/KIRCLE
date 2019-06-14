#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 21:32:44 2019

@author: galengao
"""
import sys
import warnings

import numpy as np
import pandas as pd

#import matplotlib.pyplot as plt

def matrix_difference_objective(M0, M1):
    '''Compute sum of squared differences between each element of matrices M0
    and its corresponding element in M1.'''
    return sum(sum((M0 - M1) ** 2))

def compute_belief_vector(M):
    '''Given a m-allele x n-read state matrix, compute the m-allele "belief 
    vector" of representing the state matrix's belief in each allele.'''
    V = np.sum(M, axis=1)
    return V / sum(V)

def step(M):
    '''Given a m-allele x n-read state matrix, generate the next state matrix
    using Expectation Maximization. Each read-allele belief is updated in a
    Bayesian fashion using belief vector as a prior.'''
    V = compute_belief_vector(M)
    Z = np.dot(M.T, V)
    return (M.T * V).T / Z

def run_EM(df, max_iter=1000, alpha=0.01):
    '''Run Expectation-Maximization algorithm until convergence is achieved  at
    level alpha or a maximum of n steps. Returns final belief vector after
    convergence and list of allele names it maps to'''
    # generate and normalize the initial matrix
    df.loc[:,'x'] = 1
    df0 = pd.pivot_table(df, values='x', index='sseqid', columns='qseqid', \
                         fill_value=0)
    M = np.array(df0)
    M = M / sum(M)
        
    # iteratively step our state matrix forward until convergence is acheived
    for i in range(max_iter):
        N = step(M)
        
        # if we achieve convergence before our step-limit, exit
        if matrix_difference_objective(M, N) < alpha:
            break
        
        # update state matrix
        M = N
        
    if i == max_iter-1:
        w = 'Warning: Convergence not attained after '+str(max_iter)+' steps!'
        warnings.warn(w)
        
    return compute_belief_vector(N), df0.index

if __name__ == "__main__":
#    fname = '../Workspace/KIR2DL2_blast.csv'
#    fname = 'KIR2DS2_GenotyperResults/KIR2DS2*0010101_KIR2DS2*0010108_simreads_42_blast.csv'
    fname = sys.argv[1]
    
    # load in the file to run
    df = pd.read_csv(fname, header=-1)
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',\
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = df[df.pident == 100]
    
    # run expectation maximization
    N, alleles = run_EM(df, max_iter=1000, alpha=0.0000001)
    
    # write results to file
#    outname = '_'.join(fname.split('_')[:-1]) + '_EMout.csv'
#    pd.DataFrame(N, index=alleles, columns=['p']).to_csv(outname)
    
    
#    # plot EM algorithm output in graphically appealing form
#    import matplotlib.pyplot as plt
#    plt.plot(N, 'bo')
#    plt.xticks(range(len(alleles)), alleles, rotation=30, ha='right')
#    plt.xlabel('Variant')
#    plt.ylabel('Posterior Probability')
#    plt.title('TCGA-OR-A5J2-10A; KIR3DL2', fontsize=20)
#    plt.title('KIR2DS2_Heterozygous_Simulation', fontsize=20)
