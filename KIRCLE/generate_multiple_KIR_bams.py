#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 00:19:01 2019

@author: galengao

Simple script that extracts all reads mapping to individual KIR genes.
"""
import sys
import os

import pysam

import pandas as pd

def handle_overlapping_intervals(df):
    '''Given dataframe of genomic intervals for a single KIR gene, combine
    overlapping intervals to fix double-counting of reads.'''
    def combine_overlaps(temp_tuple):
        '''Helper function to sort and then merge any overlapping intervals.'''
        temp_tuple.sort(key=lambda interval: interval[0])
        merged = [temp_tuple[0]]
        for current in temp_tuple:
            previous = merged[-1]
            if current[0] <= previous[1]:
                previous[1] = max(previous[1], current[1])
            else:
                merged.append(current)
        return merged

    # convert dataframe of intervals to list of intervals
    dfs = []
    for contig in df.chrom.unique():
        df_x = df[df.chrom == contig]
        if len(df_x) == 1: # if only 1 interval on contig, no need to combine
            dfs.append(df_x)
        else: # if multiple intervals exist that may need combining
            # extract intervals as a list of tuples
            tup = df_x[['txStart', 'txEnd']].values.tolist()
            # sort and combine tuples
            mtup = combine_overlaps(tup)
            # reformat list of tuples as a dataframe
            mtup = [[contig]+x for x in mtup[:]]
            idx = [df.index[0] for i in range(len(mtup))]
            df_out = pd.DataFrame(mtup, columns=df_x.columns, index=idx)

            dfs.append(df_out)

    return pd.concat(dfs)

def split_master_bam(infile, outname, refLocations, hg='hg19'):
    '''Given input bamfile, breaks it into multiple smaller bam files with only
    reads mapping to each KIR gene. Specify if bam is hg19 vs hg38.'''
    # read in master bam/cram/sam infile
    if infile[-4:] == '.bam':
        b = pysam.AlignmentFile(infile, "rb")
    elif infile[-4:] == '.sam':
        b = pysam.AlignmentFile(infile, "r")
    elif infile[-5:] == '.cram':
        b = pysam.AlignmentFile(infile, 'rc')
    print(infile+' loaded in.')

    # load reference table of KIR gene locations in hg19 or hg38
    df_kir = pd.read_csv(refLocations, sep='\t', index_col=0)
    # truncate table to handle overlapping regions in each KIR gene
    dfs = []
    for idx in df_kir.index.unique():
        dfs.append(handle_overlapping_intervals(df_kir[df_kir.index == idx]))
    df_kir = pd.concat(dfs).sort_index()

    # begin selecting reads in regions of interest
    newBAMfiles = {}
    for i, r in df_kir.iterrows():
        c, start, end = r['chrom'], r['txStart'], r['txEnd']

        # hg19 -- handle different alignment contig names
        if hg == 'hg19':
            if '19' not in b.references:
                if c == '19':
                    c = 'chr19'
                elif c == 'GL000209.1':
                    c = 'chr19_gl000209_random'
        # hg38 -- drop preceding 'chr' from contig names
        elif hg == 'hg38':
            c = c[3:]

        # write region's reads to gene-level miniBAM file
        fname =  outname + '_' + i + '.bam'
        # create new miniBAM if one does not exist; otherwise use existing
        if fname not in newBAMfiles:
            kirReads = pysam.AlignmentFile(fname, "wb", template=b)
            newBAMfiles[fname] = kirReads
            print(fname+' created. Writing reads...')
        else:
            kirReads = newBAMfiles[fname]

        # write reads to miniBAM file
        for read in b.fetch(c, start, end):
            kirReads.write(read)

    # finished writing. Close master bam file.
    b.close()

    # close newly made miniBAM files & sort & index the new files
    for fname in newBAMfiles.keys():
        newBAMfiles[fname].close()

        # sort and index new bam file
        pysam.sort("-o", fname, fname)
        pysam.index(fname)
        print(fname + ' sorted and indexed.')

    # return filenames of newly made miniBAMs
    return newBAMfiles.keys()

if __name__ == "__main__":
    infile = sys.argv[1]
    outname = sys.argv[2]
    refLocations = sys.argv[3]

    split_master_bam(infile, outname, f, hg='hg38')
