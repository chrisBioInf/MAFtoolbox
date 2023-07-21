#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 09:08:45 2023

@author: christopher
"""


import numpy as np
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from copy import deepcopy

from select_maf_sequences import similarity_score
from utility import check_positional_argument, read_maf, print_maf_alignment


def get_mean_gaps(records, length):
    gaps = [record.seq.count('-') / length for record in records]
    return np.mean(gaps)


def get_mean_id(records, ):
    pairwise_id = []
    
    for i in range(0, len(records)-1):
        for j in range(i, len(records)):
            seq_a = records[i].seq.replace('-', '')
            seq_b = records[j].seq.replace('-', '')
            pairwise_id.append(similarity_score(seq_a, seq_b))
            
    return np.mean(pairwise_id)


def slice_alignment(alignment, outfile, slide=40, length=120, min_id=0.25, min_seqs=2, max_gaps=0.5, max_length=120, min_length=60):
    ref_seq = str(alignment[0].seq)
    n_columns = len(ref_seq)
    result_windows = []
    
    if (n_columns <= max_length) or (n_columns <= length):
        if outfile == "":
            print_maf_alignment(alignment)
        else:
            return [alignment]
    
    slice_index = -slide
    
    while slice_index < n_columns:
        slice_index += slide
        slice_end = slice_index + length
        
        if (n_columns - slice_end) < min_length:
            break
        
        window_records = []
        
        for record in alignment:
            record_ = deepcopy(record)
            record_.seq = Seq(record.seq[slice_index:slice_end])
            record_.annotations["start"] = record.annotations["start"] + slice_index 
            record_.annotations["size"] = length -record_.seq.count("-") 
            
            if record_.annotations["size"] == 0:
                continue
            
            window_records.append(record_)
        
        if (len(window_records) < min_seqs) or (get_mean_gaps(window_records, length) > max_gaps) or (get_mean_id(window_records) < min_id):
            continue
        
        if outfile == "":
            print_maf_alignment(MultipleSeqAlignment(window_records))
        else:
            result_windows.append(MultipleSeqAlignment(window_records))
        

def window(parser):
    parser.add_option("-o","--output",action="store",type="string", default="", dest="out_file",help="MAF file to write to. If empty, results alignments are redirected to stdout.")
    parser.add_option("-s","--slide",action="store",type="int", default=40, dest="slide",help="Length of each slice step (Default: 40).")
    parser.add_option("-l","--length",action="store",type="int", default=120, dest="length",help="Length of each window (Default: 120).")
    parser.add_option("-m","--min-id",action="store",type="int", default=0.25, dest="min_id",help="Minimal mean pairwise identity (without gaps) for keeping resulting alignments, others will be dropped (Default: 0.25).")
    parser.add_option("-e","--min-seqs",action="store",type="int", default=2, dest="min_seqs",help="Minimal number of sequences in resulting alignments, those with fewer will be dropped (Default: 2).")
    parser.add_option("-g", "--max-gaps", action="store", default=0.5, type="float", dest="max_gaps", help="Maximum fraction of gaps in resulting alignments (Default: 0.5).")
    parser.add_option("-a","--max-length",action="store",type="int", default=120, dest="max_length",help="Slice only alignments longer than this value, others will be left as they are (Default: 120).")
    parser.add_option("-i","--min-length",action="store",type="int", default=120, dest="min_length",help="Leftover windows shorter than this will be discarded (Default: 60).")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    alignment_handle = read_maf(handle_)
    
    for alignment in alignment_handle:
        slice_alignment(alignment, 
                        outfile=options.out_file,
                        slide=options.slide,
                        length=options.length,
                        min_id=options.min_id,
                        min_seqs=options.min_seqs,
                        max_gaps=options.max_gaps,
                        max_length=options.max_length,
                        min_length=options.min_length
                        )
