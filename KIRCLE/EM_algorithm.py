#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 21:32:44 2019

@author: galengao
"""
from __future__ import division

import warnings

import numpy as np
import pandas as pd

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

def run_EM(df, maxIter=1000, alpha=1e-5):
    '''Run Expectation-Maximization algorithm until convergence is achieved  at
    level alpha or a maximum of n steps. Returns final belief vector after
    convergence and list of allele names it maps to'''
    # generate and normalize the initial matrix
    df.loc[:,'x'] = 1
    df0 = pd.pivot_table(df, values='x', index='sseqid', columns='qseqid', \
                         fill_value=0)
    M = np.array(df0, dtype='f')
    M = M / sum(M)
        
    # iteratively step our state matrix forward until convergence is acheived
    for i in range(maxIter):
        N = step(M)
        
        # if we achieve convergence before our step-limit, exit
        if matrix_difference_objective(M, N) < alpha:
            break
        
        # update state matrix
        M = N
        
    if i == maxIter-1:
        w = 'Warning: Convergence not attained after '+str(maxIter)+' steps!'
        warnings.warn(w)
        
    return compute_belief_vector(N), df0.index