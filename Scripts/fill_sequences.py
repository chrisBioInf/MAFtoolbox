#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 11:21:53 2023

@author: christopher
"""


from copy import deepcopy
import sys
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

from utility import read_maf, read_fasta, check_positional_argument, print_maf_alignment


def parse_mapping(map_file):
    map_ = {}
    with open(map_file, 'r') as f:
        for line in f.readlines():
            lin = line.strip("{").strip("}").strip("\n")
            if not ":" in lin:
                continue
            name, genome = lin.split(":")
            map_[name.strip()] = read_fasta(genome.strip())
            
    return map_


def fill_for_record(record, name, handle):
    if not handle:
        return record.seq
    
    seq_record = next(handle)
    
    while seq_record:
        if seq_record.id.split(".")[0] == name:
            break
        seq_record = next(handle)
        
    seq_ = str(record.seq)
    start = record.annotations.get("start")
    
    for i in range(0, len(seq_)):
        if not seq_[i] == "N":
            continue
        seq_[i] = seq_record.seq[start + i]
    
    return Seq(seq_)
    
    
def seqfill(parser):
    parser.add_option("-o", "--output", action="store", type="string", dest="output", default="", help="MAF file to write to. If empty, results alignments are redirected to stdout.")
    parser.add_option("-i", "--identity-threshold", action="store", type="float", default=1, dest="identity", help="If value is larger than 0: Check, if the not-unknown (not N) rest of the sequence fits between genome and aligned sequence with at least this fraction of its length. Serves as a sanity check if identity of genome in file and aligned genomes is not guaranteed to be identical (Default: 1)")
    parser.add_option("-m", "--mapping", action="store", type="string", dest="mapping", help="You can provide a file that maps any number of aligned sequences to a Fasta file containing the corresponding genome (Example: 'Bombus terrestris : /home/Bombus_genome.fa '), one per line.")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    required = ["mapping", ]
    
    for r in required:
        if options.__dict__[r] == None:
            print("You must pass a -%s argument." % r)
            sys.exit()
    
    genome_map = parse_mapping(options.mapping)
    alignment_handle = read_maf(handle_)
    
    for alignment in alignment_handle:
        alignment_ = []
        for record in alignment:
            if 'N' not in record.seq:
                alignment_.append(record)
                continue
            name = record.id.split(".")[0]
            record.seq = fill_for_record(record, name, genome_map.get(name, None))
            alignment_.append(deepcopy(record))
        
        if options.output == "":
            print_maf_alignment(MultipleSeqAlignment(alignment_))
            
    