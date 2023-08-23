#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:04:01 2023

@author: christopher
"""


import sys
from copy import deepcopy
import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq

from utility import read_maf, write_maf, alignments_to_stdout, max_gap_seqs, check_positional_argument


ignore_masks = True
nucleotides = ['A', 'T', 'G', 'C', '-']


def get_alignment_matrix(alignment, ignore_masks=False):
    if ignore_masks:
        column_matrix = np.array([list(record.seq.upper()) for record in alignment]).T
    else:
        column_matrix = np.array([list(record.seq) for record in alignment]).T
    return column_matrix


def slice_alignments(alignment, index_pairs, min_length=40, max_gaps=0.9):
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
            
            current_record.annotations["start"] =  record.annotations["start"] + length_previous_seq
            current_record.annotations["size"] = length_subseq
            records_.append(current_record)
        
        if len(records_) < 2:
            continue
        
        records_ = max_gap_seqs(records_, max_gaps=max_gaps, reference=True)
        if len(records_) < 2:
            continue
        sub_alignment = MultipleSeqAlignment(records_)
        result_alignments.append(sub_alignment)
    
    return result_alignments
        

def screen_for_gaps(matrix, min_gap_length=10, ignore_masks=False):
    i = 0
    j = 0
    rows_without_ref = len(matrix[0]) -1
    gap_length = 0
    index_pairs_to_cut = []
    gap_matrix = (matrix == '-')
    
    while i < len(matrix)-1:
        if (matrix[i][0] == "-"): 
            i += 1
            continue
        if not (gap_matrix[i][1:].sum() == rows_without_ref):
            i += 1
            continue
        j = i +1
        gap_length = 1
        
        while j < len(matrix):
            if (matrix[j][0] != '-') and (gap_matrix[j][1:].sum() == rows_without_ref): 
                j += 1
                gap_length += 1
            else:
                if (gap_length >= min_gap_length):
                    index_pairs_to_cut.append( (i, j) )
                gap_length = 0
                i = j +1
                break

        if gap_length >= min_gap_length:
            index_pairs_to_cut.append( (i, j) )
        if j == len(matrix):
            break
        i += 1

    return index_pairs_to_cut

        
def mask_repeat_regions(parser):
    parser.add_option("-o","--output",action="store",type="string", default="", dest="out_file",help="MAF file to write to. If empty, results alignments are redirected to stdout.")
    parser.add_option("-l","--min-length",action="store",type="int", default=40, dest="length",help="Minimal MAF block length for reporting. Shorter blocks will be dropped (Default: 40).")
    parser.add_option("-m","--min-gap-length",action="store",type="int", default=10, dest="gap_length",help="Continous gap columns shorter than this value will be ignored (Default: 10).")
    parser.add_option("-g", "--max-gaps", action="store", default=0.9, type="float", dest="max_gaps", help="All sequences with a larger gap fraction than this value will be dropped (Default: 0.9).")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    alignment_handle = read_maf(handle_)
    
    min_length = options.length
    min_gap_length = options.gap_length
    max_gap_fraction = options.max_gaps
    alignments_out = []
    
    for alignment in alignment_handle:
        if len(alignment) < 2:
            continue
        matrix = get_alignment_matrix(alignment)
        index_pairs = screen_for_gaps(matrix, min_gap_length=min_gap_length, ignore_masks=ignore_masks)
        results = slice_alignments(alignment, index_pairs, min_length=min_length, max_gaps=max_gap_fraction)
        if options.out_file != "":
            alignments_out += results
        else:
            alignments_to_stdout(results)
    
    if options.out_file != "":
        write_maf(alignments_out, options.out_file)
        
    
