#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:04:01 2023

@author: christopher
"""


import sys
from copy import deepcopy
import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq


min_length = 40
min_gap_length = 5
ignore_masks = True
max_gap_fraction = 0.75

nucleotides = ['A', 'T', 'G', 'C', '-']


def read_maf(filename):
    handle = AlignIO.parse(open(filename, 'r'), format='maf')
    return handle


def write_maf(records, filename):
    AlignIO.write(records, open(filename, 'w'), format='maf')


def print_alignment_simple(alignment):
    ref_record = alignment[0]
    print("\nAlignment:  %s    %s    %s " % (ref_record.id, ref_record.annotations["start"], ref_record.annotations["size"]))
    for record in alignment:
        print("%s: \t %s" % (record.id, record.seq))


def get_alignment_matrix(alignment, ignore_masks=False):
    if ignore_masks:
        column_matrix = np.array([list(record.seq.upper()) for record in alignment]).T
    else:
        column_matrix = np.array([list(record.seq) for record in alignment]).T
    return column_matrix


def slice_alignments(alignment, index_pairs, min_length=50):
    length = len(alignment[0].seq)
    index_pairs_to_keep_a = [0]
    index_pairs_to_keep_b = []
     
    for pair in index_pairs:
        a, b = pair
        index_pairs_to_keep_b.append(a)
        index_pairs_to_keep_a.append(b)
    
    index_pairs_to_keep_b.append(length)
    index_pairs_to_keep = list(zip(index_pairs_to_keep_a, index_pairs_to_keep_b))
    
    result_alignments = []
    print(index_pairs_to_keep)
    
    for (a, b) in index_pairs_to_keep:
        records_ = []
        if (b -a) < min_length:
            continue
        for record in alignment:
            current_record = deepcopy(record)
            current_record.seq = Seq(str(record.seq)[a:b])
            seq_ungapped = current_record.seq.replace('-', '')
            
            length_subseq = len(seq_ungapped)
            length_previous_seq = len(str(record.seq)[:a].replace('-', ''))
            
            if str(current_record.seq).count("-") / (b -a) > max_gap_fraction:
                continue
            current_record.annotations["start"] =  record.annotations["start"] + length_previous_seq
            current_record.annotations["size"] = length_subseq
            records_.append(current_record)
        
        if len(records_) < 2:
            continue
        
        sub_alignment = MultipleSeqAlignment(records_)
        print_alignment_simple(sub_alignment)
        result_alignments.append(sub_alignment)
    
    return result_alignments
        

def screen_for_gaps(matrix, min_gap_length=5, ignore_masks=False):
    i = 0
    j = 0
    rows_without_ref = len(matrix[0]) -1
    gap_length = 0
    index_pairs_to_cut = []
    
    while i < len(matrix)-1:
        if (matrix[i][0] == "-") or not (matrix[i][0].islower()):
            i += 1
            continue
        if not sum([1 for x in matrix[i][1:] if (x == "-")]) == rows_without_ref:
            i += 1
            continue
        j = i +1
        gap_length = 1
        
        while j < len(matrix):
            if (matrix[j][0] == "-") or  not (matrix[j][0].islower()):
                if gap_length >= min_gap_length:
                    index_pairs_to_cut.append( (i, j) )
                i = j +1
                break
            if sum([1 for x in matrix[j][1:] if (x == "-")]) == rows_without_ref:
                j += 1
                gap_length += 1
            else:
                if gap_length >= min_gap_length:
                    index_pairs_to_cut.append( (i, j) )
                i = j +1
                break
        i += 1

    return index_pairs_to_cut
    
    
input_ = sys.argv[1]
output_ = sys.argv[2]
handle = read_maf(input_)
alignments_out = []

for alignment in handle:
    matrix = get_alignment_matrix(alignment)
    index_pairs = screen_for_gaps(matrix, min_gap_length=min_gap_length, ignore_masks=ignore_masks)
    results = slice_alignments(alignment, index_pairs, min_length=min_length)
    alignments_out += results

write_maf(alignments_out, output_)
    
    