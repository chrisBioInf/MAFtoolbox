#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:10:16 2023

@author: christopher
"""


import pandas as pd
from Bio import AlignIO


bed_columns = ["sequence", "start", "end", "score", "strand", ]
gtf_columns = ["sequence", "source", "type", "start", "end", "score", "empty", "strand", "attribute"]


strand_dict = {
    1: "+",
    -1: "-",
    "+": "+",
    "-": "-",
    }


def read_maf(filename):
    handle = AlignIO.parse(open(filename, 'r'), format='maf')
    return handle


def write_maf(records, filename):
    AlignIO.write(records, open(filename, 'w'), format='maf')


def load_bed(filename):
    df = pd.read_csv(filename, sep="\t", header=None, names=bed_columns)
    return df


def load_gtf(filename):
    df = pd.read_csv(filename, sep="\t", header=None, names=gtf_columns)
    return df


def sortRecords(records):
    n = len(records)
    if n < 3:
        return
    swapped = False
    
    for i in range(1, n-1):
        for j in range(1, n-i-1):
            if records[j].id > records[j + 1].id:
                swapped = True
                records[j], records[j + 1] = records[j + 1], records[j]
         
        if not swapped:
            return


def print_maf_alignment(alignment):
    print("\na")
    for record in alignment:
        print("s %s \t%s \t%s \t%s \t%s \t%s" % (record.id, 
                                            record.annotations["start"], 
                                            record.annotations["size"], 
                                            strand_dict.get(record.annotations["strand"]), 
                                            record.annotations["srcSize"],
                                            record.seq))
        
def alignments_to_stdout(alignments):
    for alignment in alignments:
        print_maf_alignment(alignment)
