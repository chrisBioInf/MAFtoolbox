#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:10:16 2023

@author: christopher
"""


import sys
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq


bed_columns = ["sequence", "start", "end", "score", "strand", ]
gtf_columns = ["sequence", "source", "type", "start", "end", "score", "empty", "strand", "attribute"]


strand_dict = {
    1: "+",
    -1: "-",
    "+": "+",
    "-": "-",
    }


def read_maf(handle):
    handle = AlignIO.parse(handle, format='maf')
    return handle


def write_maf(records, filename):
    alignment = [record for record in records if record]
    
    if len(alignment) == 0:
        print("Result alignment is empty! No file written.")
    else:
        AlignIO.write(alignment, open(filename, 'w'), format='maf')


def load_bed(filename):
    df = pd.read_csv(filename, sep="\t", header=None, names=bed_columns)
    return df


def load_gtf(filename):
    df = pd.read_csv(filename, sep="\t", header=None, names=gtf_columns)
    return df


def check_positional_argument(args):
    if len(args) >= 1:
        return open(args[-1], 'r')
    else:
        return sys.stdin


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
        
        
def eliminate_consensus_gaps(records):
    ungapped_seqs = []
    seq_matrix = np.array([list(record.seq) for record in records])
    for i in range(0, len(records)):
        seq_ = ""
        for c in range(0, len(seq_matrix[i])):
            if seq_matrix[i][c] != "-":
                seq_ = seq_ + seq_matrix[i][c]
                continue
            for j in range(0, len(seq_matrix)):
                if seq_matrix[j][c] != "-":
                    seq_ = seq_ + seq_matrix[i][c]
                    break
        ungapped_seqs.append(seq_)
    
    for i in range(0, len(records)):
        records[i].seq = Seq(ungapped_seqs[i])
    
    return records


def max_gap_seqs(records, max_gaps=0, reference=True):
    start_index = 0
    if reference == True:
        start_index = 1
    length = len(records)
    
    if start_index >= length:
        return 
    
    columns = len(str(records[0].seq))
    drop_indices = set()
    
    for i in range(start_index, length):
        gaps = records[i].seq.count("-")
        gap_fraction = gaps / columns
        if (gap_fraction > max_gaps) or (gaps == columns):
            drop_indices.add(i)

    return [records[n] for n in range(0, length) if n not in drop_indices]


def print_maf_alignment(alignment):
    if not alignment:
        return
    try:
        print("\na")
        for record in alignment:
            print("s %s \t%s \t%s \t%s \t%s \t%s" % (record.id, 
                                                     record.annotations["start"], 
                                                     record.annotations["size"], 
                                                     strand_dict.get(record.annotations["strand"]), 
                                                     record.annotations["srcSize"],
                                                     record.seq))
    except BrokenPipeError:
        print("Broken pipe. Maybe check output of previous function?")
        
        
def alignments_to_stdout(alignments):
    for alignment in alignments:
        print_maf_alignment(alignment)
