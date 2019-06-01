#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 10:02:15 2019

@author: galengao
"""
import glob
import itertools

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from sklearn.metrics import confusion_matrix

import seaborn as sns

def return_sensitivity_and_specificity(dfs, i):
    tps, fps, fns, tns = {}, {}, {}, {}
    dfs_x = []
    for df, s in zip(dfs, list(sols)):
        df_x = df[[df.columns[x] for x in (2*i, 2*i+1)]]
#        df_x = df_x[df_x[df_x.columns[0]] != 'No Solution']
        truesol = tdict[s[0]]+'_'+tdict[s[1]]
        
        df_x['true'] = [truesol+'_'+str(i) for i in range(len(df_x))]
        df_x['pred'] = df_x[df_x.columns[0]] + '_' + df_x[df_x.columns[1]]
        df_x['x'] = 1
        df_x = df_x.pivot_table(index='true', columns='pred', values='x', \
                                fill_value=0)
        dfs_x.append(df_x)
        
        if truesol in df_x.columns:
            tps[truesol] = sum(df_x[truesol])
            fns[truesol] = len(df_x) - sum(df_x[truesol])
        else:
            tps[truesol] = 0
            fns[truesol] = len(df_x)
        
    df0 = pd.concat(dfs_x).fillna(value=0)
    for s in sols:
        truesol = tdict[s[0]]+'_'+tdict[s[1]]
        if truesol in df0.columns:
            fps[truesol] = sum(df0[truesol]) - tps[truesol]
            tns[truesol] = len(df0) - sum(df0[truesol]) - fns[truesol]
        else:
            fps[truesol] = 0
            tns[truesol] = len(df0) - fns[truesol]
    
    tp, fn, fp, tn = sum(tps.values()), sum(fns.values()), sum(fps.values()), sum(tns.values())
    
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    
    return sens, spec

def plot_confusion_matrix(y_true, y_pred, classes, normalize=False, title=None, cmap=plt.cm.Blues, text=True):
    """Print & plot confusion matrix. Normalize by setting `normalize=True`."""
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
#    classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=50, ha="right",
             rotation_mode="anchor", fontsize=7)
    plt.setp(ax.get_yticklabels(), fontsize=7)

    # Loop over data dimensions and create text annotations.
    if text:
        fmt = '.2f' if normalize else 'd'
        thresh = cm.max() / 2.
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, format(cm[i, j], fmt),
                        ha="center", va="center", fontsize=8,
                        color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax
    
def kir_genotype_confusion_matrix(dfs, alleles, i, text=True):
    '''Use simulation data processed output to generate confusion matrix. Use
    specified column 2i and 2i+1.'''
    sols = list(itertools.combinations_with_replacement(alleles, r=2))

    transdict = {tdict[s[0]]+'_'+tdict[s[1]]:j for j,s in enumerate(sols)}
    for j, s in enumerate(sols):
        transdict[tdict[s[1]]+'_'+tdict[s[0]]] = j
    transdict['No Solution_No Solution'] = -1
#    print('KIR:KIR00429_KIR:KIR00636' in transdict)
#    print('KIR:KIR00636_KIR:KIR00429' in transdict)
    trueGenotypes, predGenotypes = [], []
    for s, df in zip(sols, dfs):
        solName = tdict[s[0]] + '_' + tdict[s[1]]
        solNameSeries = df[df.columns[2*i]] + '_' + df[df.columns[2*i+1]]
        print(solName)
        print(solNameSeries[0])
        predGenotypes.append([transdict[s] for s in solNameSeries])
        trueGenotypes.append([transdict[solName] for i in range(len(df))])
    
    predGenotypes = np.ravel(predGenotypes)
    trueGenotypes = np.ravel(trueGenotypes)
    
    
    np.set_printoptions(precision=2)
    
    # Plot non-normalized confusion matrix
    sols = list(itertools.combinations_with_replacement(alleles, r=2))
    class_names  = ['No Solution_No Solution'] + [tdict[s[0]]+'_'+tdict[s[1]] for s in sols]
    plot_confusion_matrix(trueGenotypes, predGenotypes, classes=class_names,
                          title='Confusion matrix, without normalization', text=text)
    plt.gcf().set_size_inches((20,15))
    plt.savefig('ConfusionMatrix_t='+str((i+1)*0.02)+'.png', bbox_inches='tight')
#    plt.clf()
    
def ROC_curve(dfs):
    '''Plot ROC curve for data in dataframes.'''
    # rearrange dataframes
    ys, xs = [], []
    for i in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]:
#    for i in [0,1,2,3,4,5,6,7,8,9]:
        print(i)
        sens, spec = return_sensitivity_and_specificity(dfs, i)
        ys.append(sens)
        xs.append(1-spec)
    
    plt.plot(xs, ys, 'rx')
    #plt.plot([min(xs), max(xs)], [min(xs), max(xs)], 'k-')
    plt.xlabel('1 - Specificity')
    plt.ylabel('Sensitivity')

    
# tabulate all potential 2-allele solutions
kgene = 'KIR2DS3'
df_kirvars = pd.read_csv('../KIR_variants.txt', sep='\t')
alleles = df_kirvars[kgene].dropna()
sols = list(itertools.combinations_with_replacement(alleles, r=2))

# KIR2DS1
#tdict = {'KIR2DS1*0020101':'KIR:KIR00034', 'KIR2DS1*0020103':'KIR:KIR00296',\
#         'KIR2DS1*0020102':'KIR:KIR00298', 'KIR2DS1*0020104':'KIR:KIR00426',\
#         'KIR2DS1*0020105':'KIR:KIR00427', 'KIR2DS1*0020106':'KIR:KIR00428'}

# KIR2DS2
#tdict = {'KIR2DS2*0010101':'KIR:KIR00037', 'KIR2DS2*0010102':'KIR:KIR00302',\
#         'KIR2DS2*0010103':'KIR:KIR00270', 'KIR2DS2*0010104':'KIR:KIR00336',\
#         'KIR2DS2*0010105':'KIR:KIR00446', 'KIR2DS2*0010106':'KIR:KIR00447',\
#         'KIR2DS2*0010107':'KIR:KIR00448', 'KIR2DS2*0010108':'KIR:KIR00449',\
#         'KIR2DS2*0010109':'KIR:KIR00450', 'KIR2DS2*0010110':'KIR:KIR00451',\
#         'KIR2DS2*0010111':'KIR:KIR00623', 'KIR2DS2*0010112':'KIR:KIR00635'}

# KIR2DS5
#tdict = {'KIR2DS5*0020101': 'KIR:KIR00050', 'KIR2DS5*0020102':'KIR:KIR00299',\
#         'KIR2DS5*0020103': 'KIR:KIR00297', 'KIR2DS5*0020104':'KIR:KIR00433',\
#         'KIR2DS5*006':'KIR:KIR00276', 'KIR2DS5*007':'KIR:KIR00277',\
#         'KIR2DS5*010':'KIR:KIR00432'}

# KIR2DS3
tdict = {'KIR2DS3*0010301': 'KIR:KIR00044', 'KIR2DS3*0010302': 'KIR:KIR00636',\
         'KIR2DS3*0020101': 'KIR:KIR00131', 'KIR2DS3*0020102': 'KIR:KIR00429',\
         'KIR2DS3*0020103': 'KIR:KIR00652'}

# KIR2DL5A
#tdict = {'KIR2DL5A*0010101':'KIR:KIR00029', 'KIR2DL5A*0010102':'KIR:KIR00300',\
#         'KIR2DL5A*0010103':'KIR:KIR00439', 'KIR2DL5A*00102':'KIR:KIR00259',\
#         'KIR2DL5A*00103':'KIR:KIR00260', 'KIR2DL5A*00104':'KIR:KIR00262',\
#         'KIR2DL5A*00105':'KIR:KIR00263', 'KIR2DL5A*0050101':'KIR:KIR00095',\
#         'KIR2DL5A*0050102':'KIR:KIR00249', 'KIR2DL5A*0050103':'KIR:KIR00438',\
#         'KIR2DL5A*0050104':'KIR:KIR00440', 'KIR2DL5A*01201':'KIR:KIR00391',\
#         'KIR2DL5A*01202':'KIR:KIR00389'}

# KIR2DL2
#tdict = {'KIR2DL2*0010101':'KIR:KIR00010', 'KIR2DL2*0010102':'KIR:KIR00405',\
#         'KIR2DL2*0010103':'KIR:KIR00407', 'KIR2DL2*0010104':'KIR:KIR00410',\
#         'KIR2DL2*0010105':'KIR:KIR00618', 'KIR2DL2*0010106':'KIR:KIR00645',\
#         'KIR2DL2*0010107':'KIR:KIR00630', 'KIR2DL2*0030101':'KIR:KIR00012',\
#         'KIR2DL2*0030102':'KIR:KIR00316', 'KIR2DL2*0030103':'KIR:KIR00335',\
#         'KIR2DL2*0030104':'KIR:KIR00406', 'KIR2DL2*0030105':'KIR:KIR00408',\
#         'KIR2DL2*0030106':'KIR:KIR00409', 'KIR2DL2*0030107':'KIR:KIR00411'}

# KIR3DS1
#tdict = {'KIR3DS1*0130101':'KIR:KIR00085', 'KIR3DS1*0130102':'KIR:KIR00489',\
#         'KIR3DS1*0130103':'KIR:KIR00488', 'KIR3DS1*0130104':'KIR:KIR00708',\
#         'KIR3DS1*0130105':'KIR:KIR00709', 'KIR3DS1*0130106':'KIR:KIR00710',\
#         'KIR3DS1*0130107':'KIR:KIR00711', 'KIR3DS1*0130108':'KIR:KIR00713',\
#         'KIR3DS1*0130109':'KIR:KIR00715','KIR3DS1*0130110':'KIR:KIR00717',\
#         'KIR3DS1*0130111':'KIR:KIR00719', 'KIR3DS1*014':'KIR:KIR00086',\
#         'KIR3DS1*055':'KIR:KIR00243','KIR3DS1*078':'KIR:KIR00751'}

# load solutions
#stem = 'KIR2DS1_GenotyperResults/'
stem = 'runs_'+kgene+'_bsp70/'
tail = 'Genotyper_results_different_thresholds.csv'
dfs = [pd.read_csv(f, index_col=0) for f in [stem+s[0]+'_'+s[1]+tail for s in sols]]
#df = pd.read_csv('KIR2DS1*0020101_KIR2DS1*0020101Genotyper_results_different_thresholds.csv', index_col=0)


# generate ROC curve
#ROC_curve(dfs)


# generate confusion matrix
#kir_genotype_confusion_matrix(dfs, alleles, 0, text=True)


# visualize distance matrix
data = []
with open('../clustal_omega_results/percent_identity_matrices/KIR2DL3_pim.txt', 'r') as f:
    for line in f:
        if line[0] != '#':
            data.append(line.split())

ticks = np.array(data[1:]).T[1].T
sns.heatmap(np.array(data[1:]).T[2:].T.astype(np.float), vmin=99.95, cmap='Reds',\
            xticklabels=ticks, yticklabels=ticks)




#tps, fps, fns, tns = 0, 0, 0, 0
#for s in sols:
#    a1_truth = np.array(df[a1] == transdict[s[0]])
#    a2_truth = np.array(df[a2] == transdict[s[1]])
#    sens = sum( np.logical_and(a1_truth, a2_truth) )
#    tps += sens
#
#ts = np.linspace(0.02, 0.2, 10)
#
#
#sensitivities, specificities = {}, {}
#for a1, a2 in zip(*[iter(df.columns)]*2):
#    a1_truth = np.array(df[a1] == 'KIR:KIR00034')
#    a2_truth = np.array(df[a2] == 'KIR:KIR00034')
#    sens = sum( np.logical_and(a1_truth, a2_truth) ) / len(df)
#    
#    print(sens)