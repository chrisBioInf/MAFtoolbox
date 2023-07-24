#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:31:30 2023

@author: christopher
"""


from copy import deepcopy
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
import sys

from utility import write_maf, alignments_to_stdout, check_positional_argument, read_maf
from highlight_regions import load_bed_with_range


def count_nucleotides(seq):
    count = 0
    for c in seq:
        if c == "-":
            continue
        count += 1
    
    return count


def get_subrecord(record, start, end):
    seq_ = str(record.seq)
    extracted_region = seq_[start:end]
    new_length = count_nucleotides(extracted_region)
    new_start = record.annotations.get("start") + count_nucleotides(seq_[:start])
    
    new_record = deepcopy(record)
    new_record.seq = Seq(extracted_region)
    new_record.annotations["start"] = new_start
    new_record.annotations["size"] = new_length
    
    return new_record


def get_subblock(alignment, df, sense, antisense, min_seqs):
    start = int(alignment[0].annotations["start"]) 
    size = int(alignment[0].annotations["size"])
    end = start + size
    start_end = set((x for x in range(start, end)))
    
    alignments = []
    
    for i in range(0, len(df)):
        if len(df["range"].iloc[i].intersection(start_end)) == 0:
            continue
        
        records = []
        
        anno_size = df["end"].iloc[i] -df["start"].iloc[i] 
        extracted_block_size = anno_size + sense + antisense
        offset_start = df["start"].iloc[i] -start -antisense
        offset_end = offset_start + extracted_block_size
        
        if offset_start < 0:
            offset_start = 0
            
        section1 = ''
        extracted_section = ''
        section2 = '' 
        n = 0
        gaps = 0
        ref_seq = str(alignment[0].seq) 
        
        while n < len(ref_seq):
            char_ = ref_seq[n]
            if (n -gaps >= offset_start):
                if (n-gaps >= offset_end):
                    section2 += char_
                else:
                    extracted_section += char_
            else:
                section1 += char_
        
            if char_ == "-":
                gaps += 1
            n += 1
        
        start_index = len(section1)
        end_index = len(section1) + len(extracted_section)
        length = end_index - start_index
        
        for record in alignment:
            record_ = get_subrecord(record, start_index, end_index)
            if record_.seq.count('-') < length:
                records.append(record_)
        
        if len(records) >= min_seqs:
            alignments.append(MultipleSeqAlignment(records))
    
    return alignments


def extract_blocks(handle_, annotation, sense=0, antisense=0, output="", subset=[], min_seqs=2):
    output_blocks = []
    alignment_handle = read_maf(handle_)
    subsections = {}
    
    for name in annotation["sequence"].unique():
        subsections[name] = annotation[annotation["sequence"] == name]
    
    for alignment in alignment_handle:
        ref_id = str(alignment[0].id)
        annotations_ = subsections.get(ref_id, [])
        
        if len(annotations_) == 0:
            continue
        
        if subset:
            subset.add(ref_id.split(".")[0].lower())
            alignment_ = [record for record in alignment if record.id.split(".")[0].lower() in subset]
            extracted_alignments = get_subblock(alignment_, annotations_, sense, antisense, min_seqs)
        else:
            extracted_alignments = get_subblock(alignment, annotations_, sense, antisense, min_seqs)
        
        if output == "":
            alignments_to_stdout(extracted_alignments)
        else:
            output_blocks += extracted_alignments
            
    return output_blocks


def extract_alignment(parser):
    parser.add_option("-b","--bed",action="store", type="string", dest="bed", help="Bed file with genomic coordinates to extract (Required).")
    parser.add_option("-u","--subset",action="store",type="string",dest="subset",default="",help="You may select sequences by entering a comma separated list of sequence names here (Example: Apis_mellifera,Bombus_terrestris). All other sequences will be discarded, but reference will always be kept.")
    parser.add_option("-s", "--sense",action="store",type="int",dest="sense",default=0,help="Add an overhang of this many nucleotides in sense (+) direction of reference strand (Default: 0).")
    parser.add_option("-n", "--antisense",action="store",type="int",dest="antisense",default=0,help="Add an overhang of this many nucleotides in antisense (-) direction of reference strand (Default: 0).")
    parser.add_option("-o","--output",action="store",type="string", default="", dest="output", help="MAF file to write to. If empty, results alignments are redirected to stdout.")
    parser.add_option("-r", "--remove-duplicates",action="store_true",default=False,dest="remove_duplicates",help="Should identical coordinates be filtered out? ()")
    parser.add_option("-m","--min-seqs",action="store",type="int", default=2, dest="min_seqs", help="Only report blocks with at least this many remaining sequences. All-gap sequences will be dropped. (Default: 2).")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    required = ["bed", ]
    
    for r in required:
        if options.__dict__[r] == None:
            print("You must pass a --%s argument." % r)
            sys.exit()
            
    annotations = load_bed_with_range(options.bed)
    
    if options.remove_duplicates == True:
        annotations.drop_duplicates(subset=["start", "end"], inplace=True)
        
    if options.subset != "":
        subset = set(x.lower() for x in options.subset.split(","))
    else:
        subset = []
    
    extracted_blocks = extract_blocks(handle_, annotations, options.sense, options.antisense, options.output, subset, options.min_seqs)
    
    if len(extracted_blocks) > 0:
        write_maf(extracted_blocks, options.output)
        
    