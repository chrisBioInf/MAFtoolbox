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
from Bio.Seq import Seq


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


def write_maf(records, filename):
    AlignIO.write(records, open(filename, 'w'), format='maf')


def print_alignment_simple(alignment):
    ref_record = alignment[0]
    print("\nAlignment:  %s    %s  " % (ref_record.id, ref_record.annotations["start"]))
    for record in alignment:
        print("%s: \t %s" % (record.id, record.seq))
        
        
def eliminate_consensus_gaps(records):
    ungapped_seqs = []
    for i in range(0, len(records)):
        seq = str(records[i].seq)
        seq_ = ""
        for c in range(0, len(seq)):
            if seq[c] != "-":
                seq_ = seq_ + seq[c]
                continue
            for j in range(0, len(records)):
                seq_2 = str(records[j].seq)
                if seq_2[c] != "-":
                    seq_ = seq_ + seq[c]
                    break
        ungapped_seqs.append(seq_)
    
    for i in range(0, len(records)):
        records[i].seq = Seq(ungapped_seqs[i])
    
    return records


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
        gap_fraction = consensus_seq.count('-') / length
        
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
            
        current_records = eliminate_consensus_gaps(current_records)
        
    if alignment_discarded:
        return MultipleSeqAlignment([alignment[0]])
    
    return MultipleSeqAlignment(current_records)
    

def main():
    input_ = sys.argv[1]
    output_ = sys.argv[2]
    handle = read_maf(input_)
    output_alignments = []
    
    for alignment in handle:
        alignment = filter_records_by_similarity(alignment)
        if len(alignment) > 1:
            print_alignment_simple(alignment)
            output_alignments.append(alignment)
    
    write_maf(output_alignments, output_)
    
    
main()
