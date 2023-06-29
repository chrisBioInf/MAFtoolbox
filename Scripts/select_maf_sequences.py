#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 11:36:23 2023

@author: christopher
"""


import sys
import numpy as np
from difflib import SequenceMatcher
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


ignore_masks = True
keep_reference = True
id_threshold = 0.8
min_seqs = 6
max_gaps = 0.5

nucleotides = ['A', 'T', 'G', 'C', '-']


def similarity_score(a, b):
    return SequenceMatcher(None, a, b).ratio()


def read_maf(filename):
    handle = AlignIO.parse(open(filename, 'r'), format='maf')
    return handle


def print_alignment_simple(alignment):
    ref_record = alignment[0]
    print("\nAlignment:  %s    %s  " % (ref_record.id, ref_record.annotations["start"]))
    for record in alignment:
        print("%s: \t %s" % (record.id, record.seq))
        
        
def remove_consensus_gaps(alignment):
    pass


def calculate_consensus_sequence(records):
    consensus_nucleotides = []
    
    if ignore_masks:
        seq_columns = np.array([list(record.seq.upper()) for record in records]).T
    else:
        seq_columns = np.array([list(record.seq) for record in records]).T

    for i in range(0, len(seq_columns)):
        unique, counts = np.unique(seq_columns[i], return_counts=True)
        nt_dict = dict(zip(unique, counts))
        nt_counts = [nt_dict.get('A', 0), nt_dict.get('T', 0),
                     nt_dict.get('G', 0), nt_dict.get('C', 0),
                     nt_dict.get('-', 0)]
        consensus_nucleotides.append(nucleotides[np.argmax(nt_counts)])
    
    consensus_seq = str().join(consensus_nucleotides)
    return consensus_seq


def filter_records_by_similarity(alignment):
    length = len(alignment)
    current_records = [record for record in alignment]
    length = len(current_records[0].seq)
    
    if len(alignment) < 3:
        return alignment
    
    alignment_discarded = False
    alignment_finalized = False
    
    while not alignment_finalized:
        consensus_seq = calculate_consensus_sequence(current_records)
        print(consensus_seq)
        gap_fraction = consensus_seq.count('-') / length
        print(gap_fraction)
        
        if gap_fraction > max_gaps:
            alignment_discarded = True
            break
        
        similarity_column = [similarity_score(record.seq, consensus_seq) for record in current_records]
        
        minimal_indices = sorted(np.argpartition(similarity_column, 2)[:2])
        minimal_index = minimal_indices[0]
        if keep_reference and (minimal_index == 0):
            minimal_index = minimal_indices[1]
        current_records.pop(minimal_index)
        
        if len(current_records) <= min_seqs:
            alignment_finalized = True 
            
        if np.mean(similarity_column) > id_threshold:
            alignment_finalized = True
        
    if alignment_discarded:
        return MultipleSeqAlignment([alignment[0]])
    
    return MultipleSeqAlignment(current_records)
    

def main():
    input_ = sys.argv[1]
    handle = read_maf(input_)
    
    for alignment in handle:
        alignment = filter_records_by_similarity(alignment)
        print_alignment_simple(alignment)
    
main()
