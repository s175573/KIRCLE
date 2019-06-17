#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 00:50:09 2019

@author: galengao
"""
import numpy as np
import pandas as pd

import kir_genotyper.EM_algorithm as EM
import kir_genotyper.bootstrapper as bs
import kir_genotyper.genotyper as gt

def test_EM_algorithm():
    '''Test EM algorithm with set numbers.'''
    M = np.array([[1/3, 1/2, 0],
                  [1/3, 0, 0],
                  [1/3, 1/2, 1]])
    
    M_prime = np.array([[5/18, 5/16, 0], 
                        [2/18, 0, 0],
                        [11/18, 11/16, 1]])
    
    # test EM belief vector function
    if EM.compute_belief_vector(M) == np.array([5/6, 2/6, 11/6]):
        print('EM belief vector function passed test 1/1.')
    else:
        print('EM belief vector function failed test 1/1.')
    
    # test EM step function
    if EM.step(M) == M_prime:
        print('EM step function passed test 1/1.')
    else:
        print('EM step function failed test 1/1.')
    
def test_bootstrapper():
    '''Test bootstrapper.'''
    
    # test bootstrapper booststrap_readnames function
    passed = 0
    uniqueReads = ['a,', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
    for p in [0.2, 0.4, 0.6, 0.8, 1]:
        bsreads = bs.bootstrap_readnames(uniqueReads, part=p)
        if len(bsreads) == p * 10 + 1:
            passed += 1
            
    print('Bootstrapper bootstrap_readnames function passed ' + str(passed) + \
          '/10 tests.')
        
    # test bootstrapper bootstrap function
    ######################################
    ######################################
    ######################################
    
    # test case of empty dataframe
    df = pd.DataFrame([])
    for k in kgenes:
        P, alleles = bs.bootstrap_EM(df, k)
        
        # pull alleles from file
        df_kirvars = pd.read_csv('../KIR_variants.txt', sep='\t')
        sortedAlleles = sorted([tdict[a] for a in df_kirvars[k].dropna()])
    
        passed = 0
        if np.all(alleles == sortedAlleles):
            if np.all(P == np.array([0 for x in sortedAlleles])):
                passed += 1
        
    print('Empty Dataframe passed bootstrapper tests: ' + \
          str(len(kgenes) - passed) + '/' + str(len(kgenes)))
    
    # test cases of full homozygous and heterozygous dataframes
    np.random.seed(42)
    
    answers = [[0, 0, ...],
               [0, 0, ...]]
    for l, key in zip(['homozygous','heterozygous'], answers):
        df = pd.read_csv(l + '.csv', header=-1)
        df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',\
                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    
        P, alleles = bs.bootstrap_EM(df, 'KIR3DL1', part=0.3, n_boot=5, max_iter=100)
        
        passed = False
        if np.all(alleles == sorted([tdict[a] for a in df_kirvars['KIR3DL1'].dropna()])):
            if np.all(P == np.array(key)):
                print('Full dataframe of ' + l + ' simulated genome passed.')
                

def test_kir_genotyper():
    '''Test KIR genotyper.'''
    gt.genotyper(df)
    
        
kgenes = ['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5A', 'KIR2DL5B',\
          'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5', 'KIR3DS1',\
          'KIR3DL1', 'KIR3DL2', 'KIR2DP1', 'KIR3DP1']

df_tr = pd.read_csv('../KIRallele_translation.txt', sep='\t', index_col=0)
tdict = df_tr['Colon_Name'].to_dict()

print(' --> Testing EM algorithm <-- ')
test_EM_algorithm()
print(' ---------------------------- ')
print(' ')
print(' ')

print(' --> Testing bootstrapper <-- ')
test_bootstrapper()
print(' ---------------------------- ')
print(' ')
print(' ')

print(' --> Testing KIR genotyper <-- ')
test_kir_genotyper()
print(' ---------------------------- ')