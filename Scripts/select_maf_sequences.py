#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 11:36:23 2023

@author: christopher
"""


import numpy as np
from difflib import SequenceMatcher
from Bio.Align import MultipleSeqAlignment

from utility import read_maf, write_maf, print_maf_alignment, eliminate_consensus_gaps, max_gap_seqs, check_positional_argument


ignore_masks = True

nucleotides = ['A', 'T', 'G', 'C', '-']


def similarity_score(a, b):
    return SequenceMatcher(None, a, b).ratio()


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


def filter_records_by_similarity(alignment, options):
    keep_reference = options.reference
    id_threshold = options.id_threshold
    min_seqs = options.min_seqs
    current_records = [record for record in alignment]
    
    if len(alignment) < 3:
        return alignment
    
    alignment_discarded = False
    alignment_finalized = False
    
    while not alignment_finalized:
        consensus_seq = calculate_consensus_sequence(current_records)
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


def select_seqs(parser):
    parser.add_option("-o","--output",action="store",type="string", default="", dest="out_file",help="MAF file to write to. If empty, results alignments are redirected to stdout.")
    parser.add_option("-r", "--no-reference", action="store_false", default=True, dest="reference", help="Should the first sequence always be considered the reference? (Default: True)")
    parser.add_option("-p","--id-threshold",action="store", type="float", default=0.8, dest="id_threshold", help="No further sequences are removed from alignment if average pairwise identity to the consensus sequence is equal to or larger than this value (Default: 0.8).")
    parser.add_option("-m","--min-seqs",action="store",type="int", default=6, dest="min_seqs", help="No further sequences will be removed if input alignment reaches this number or fewer sequences (Default: 6).")
    parser.add_option("-g", "--max-gaps", action="store", default=0.9, type="float", dest="max_gaps", help="All sequences with a larger gap fraction than this value will be dropped (Default: 0.9).")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    alignment_handle = read_maf(handle_)
    output_alignments = []
    
    for alignment in alignment_handle:
        alignment = max_gap_seqs(alignment, max_gaps=options.max_gaps, reference=options.reference)
        alignment = filter_records_by_similarity(alignment, options)
        if len(alignment) > 1:
            if options.out_file == "":
                print_maf_alignment(alignment)
            else:
                output_alignments.append(alignment)
    
    if len(output_alignments) > 0:
        write_maf(output_alignments, options.out_file)
    
