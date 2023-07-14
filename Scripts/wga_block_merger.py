#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:52:16 2023

@author: christopher
"""


import sys
from Bio import Seq
import numpy as np

from utility import write_maf, read_maf, print_maf_alignment, sortRecords, eliminate_consensus_gaps, max_gap_seqs, check_positional_argument


def coordinate_distance(end1, start2):
    return start2 - end1


def pairwise_sequence_identity(seq1, seq2):
    count = 0
    for i in range(0, len(seq1)):
        if seq1[i] == seq2[i]:
            count += 1
    return count / len(seq1)
 
 
def local_species_consensus(block1, block2):
    block1_ids = set([str(seq.id) for seq in block1])
    block2_ids = set([str(seq.id) for seq in block2])
    consensus_species = len(block1_ids.intersection(block2_ids)) + len(block2_ids.intersection(block1_ids)) 
    consensus_score  = consensus_species / (len(block1_ids) + len(block2_ids))
    return consensus_species, consensus_score


def concat_with_bridge(seq1, seq2, offset, max_offset):
    n_gaps = max_offset - offset
    seq_ = Seq.Seq(str(seq1) + "N"*offset + '-'*n_gaps + str(seq2)) 
    return seq_


def filter_by_seq_length(records, reference=True):
    length_dict = {str(record.id) : len(record.seq.replace('-', '')) for record in records}
    record_dict = {str(record.id) : record for record in records}
    mu = np.mean([x for x in length_dict.values()])
    std = np.std([x for x in length_dict.values()])
    
    count = 0
    
    for record in records:
        count += 1
        if (count == 1) and reference:
            continue
        
        x = str(record.id)
        if std > 0:
            z = abs(length_dict.get(x) - mu) / std
        else:
            z = 0
        if z > 2:
            record_dict.pop(x)
    
    return [record for record in record_dict.values()]


def merge_blocks(block1, block2, reference=True, offset_threshold=0):
    block1_ids = [str(record.id) for record in block1]
    block2_ids = [str(record.id) for record in block2]
    
    if len(block1_ids) >= len(block2_ids):
        consensus_blocks = set(block1_ids).intersection(block2_ids)
    else:
        consensus_blocks = set(block2_ids).intersection(block1_ids)
    
    block1_records = [record for record in block1 if record.id in consensus_blocks]
    block2_records = [record for record in block2 if record.id in consensus_blocks]
    sortRecords(block1_records)
    sortRecords(block2_records)
    
    merged_records = []
    offsets = []
    
    for i in range(0, len(block1_records)):
        offsets.append( coordinate_distance(block1_records[i].annotations["start"] + block1_records[i].annotations["size"],
                                            block2_records[i].annotations["start"] ))
    max_offset = max(offsets)
    
    for i in range(0, len(block1_records)):
        record_, record2 = block1_records[i], block2_records[i]
        strand1, strand2 = record_.annotations["strand"], record2.annotations["strand"]
        start1, end1 = record_.annotations["start"], record_.annotations["start"] + record_.annotations["size"]
        start2, end2 = record2.annotations["start"], record2.annotations["start"] + record2.annotations["size"]
        offset = offsets[i]
        if ((offset > offset_threshold) or (strand1 != strand2)):
            continue
        record_.seq = concat_with_bridge(record_.seq, record2.seq, offset, max_offset)
        record_.annotations["size"] = (end1 - start1) + offset + (end2 - start2)
        merged_records.append(record_)
    
    return merged_records


def check_block_viability(block1, block2, species_consensus_threshold, block_distance_threshold, block_length_threshold, min_seqs=2):
    merge_flag = True
    reference1 = block1[0]
    reference2 = block2[0]
    start1, end1 = reference1.annotations["start"], reference1.annotations["start"] + reference1.annotations["size"]
    start2, end2 = reference2.annotations["start"], reference2.annotations["start"] + reference2.annotations["size"]
    consensus_species, consensus_score = local_species_consensus(block1, block2)
    
    if coordinate_distance(end1, start2) > block_distance_threshold:
        merge_flag = False
    elif reference1.annotations["strand"] != reference2.annotations["strand"]:
        merge_flag = False
    elif (len(reference1.seq) + len(reference2.seq)) > block_length_threshold:
        merge_flag = False 
    elif consensus_score < species_consensus_threshold:
        merge_flag = False
    elif consensus_species < min_seqs:
        merge_flag = False
    elif reference1.id != reference2.id:
        merge_flag = False
         
    return merge_flag


def merge(parser):
    parser.add_option("-o", "--output", action="store", type="string", dest="out_file", default="", help="MAF file to write to. If empty, results alignments are redirected to stdout.")
    parser.add_option("-r", "--no-reference", action="store_false", default=True, dest="reference", help="Set this flag if the first sequence should NOT be considered as reference.")
    parser.add_option("-s", "--species-consensus", action="store", type="float", default=0.75, dest="species_consensus", help="Minimal consensus between neighboring blocks for merging (Default: 0.75).")
    parser.add_option("-d", "--max-distance", action="store", default=0, type="int", dest="distance", help="Maximum distance between genomic coordinates of sequences for merging of neighboring blocks (Default: 0).")
    parser.add_option("-l", "--max-length", action="store", default=1000, type="int", dest="length", help="Merged alignment blocks will not be extended past this block length (Default: 1000).")
    parser.add_option("-g", "--max-gaps", action="store", default=0.9, type="float", dest="max_gaps", help="All sequences with a larger gap fraction than this value will be dropped (Default: 0.9).")
    parser.add_option("-m", "--min-seqs", action="store", default=2, type="int", dest="min_seqs", help="No merging will happen, if the blocks share this few or less sequences (Default: 2).")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    alignments = read_maf(handle_)
    merged_alignments = []
    block1 = next(alignments)
    merged_blocks_n = []
    local_merges = 1
    
    while block1 != None:
        block2 = next(alignments, None)
        
        if block2 == None:
            records = eliminate_consensus_gaps(block1)
            records = max_gap_seqs(records, options.max_gaps, options.reference)
            merged_blocks_n.append(local_merges)
            if options.out_file == "":
                print_maf_alignment(records)
            else:
                merged_alignments.append(records)
            break
                
        merge_flag = check_block_viability(block1, block2, 
                                           options.species_consensus,
                                           options.distance,
                                           options.length)
        if merge_flag == True: 
            block1 = merge_blocks(block1, block2, options.reference, 
                                  offset_threshold=options.distance)
            local_merges += 1
        else:
            records = eliminate_consensus_gaps(block1)
            records = max_gap_seqs(records, options.max_gaps, options.reference)
            merged_blocks_n.append(local_merges)
            local_merges = 1
            block1 = block2
            if options.out_file == "":
                print_maf_alignment(records)
            else:
                merged_alignments.append(records)

    if options.out_file != "":
        write_maf(merged_alignments, options.out_file)
    
    
