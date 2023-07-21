#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:10:29 2023

@author: christopher
"""


import sys
import pandas as pd

from utility import strand_dict, check_positional_argument, read_maf


columns = ["sequence", "start", "end", "score", "strand", ]

colors_dict = {
    "BLUE" : '\033[94m',
    "GREEN" : '\033[92m',
    "RED" : '\033[91m',
    "YELLOW" : '\033[93m',
    "END" : '\033[0m',
}


def load_bed_with_range(filename):
    df = pd.read_csv(filename, sep="\t", header=None, names=columns)
    ranges = []
    
    for i in range(0, len(df)):
        start = df["start"].iloc[i]
        end = df["end"].iloc[i]
        ranges.append(frozenset((x for x in range(start, end))))
    
    df["range"] = ranges
    return df 


def print_highlighted_sequence(record, color_start_index, color_end_index, 
                               sense_overhang_start, sense_overhang_end, 
                               antisense_overhang_start, antisense_overhang_end, 
                               color, overhang_color):
    seq_ = str(record.seq)
    start = record.annotations.get("start")
    size = record.annotations.get("size")
    strand = strand_dict.get(record.annotations.get("strand"), "+")
    srcsize = record.annotations.get("srcSize")
    
    section1 = seq_[:antisense_overhang_start]
    overhang1 = seq_[antisense_overhang_start:antisense_overhang_end]
    highlight_region = seq_[color_start_index:color_end_index]
    overhang2 = seq_[sense_overhang_start:sense_overhang_end]
    section2 = seq_[sense_overhang_end:]

    print('s %s \t%s \t%s \t%s \t%s \t %s%s%s%s%s%s%s%s%s%s%s' % (
                                                record.id, start, size, strand, srcsize,
                                                section1, 
                                                colors_dict.get(overhang_color, '\033[92m'), 
                                                overhang1, 
                                                colors_dict.get("END"), 
                                                colors_dict.get(color, '\033[92m'), 
                                                highlight_region, 
                                                colors_dict.get("END"), 
                                                colors_dict.get(overhang_color, '\033[92m'), 
                                                overhang2, 
                                                colors_dict.get("END"), 
                                                section2)
          )


def print_highlighted_alignment(alignment, df, sense=0, antisense=0, color="GREEN", overhang_color="GREEN"):
    start = int(alignment[0].annotations["start"])
    size = int(alignment[0].annotations["size"])
    end = start + size
    start_end = set((x for x in range(start, end)))
    
    for i in range(0, len(df)):
        if len(df["range"].iloc[i].intersection(start_end)) == 0:
            continue
        anno_size = df["end"].iloc[i] -df["start"].iloc[i]
        offset_start = df["start"].iloc[i] -start
        offset_end = offset_start + anno_size
        
        if offset_start < 0:
            offset_start = 0
        
        section1 = ''
        overhang1 = ''
        highlight_section = ''
        overhang2 = ''
        section2 = '' 
        n = 0
        gaps = 0
        ref_seq = str(alignment[0].seq) 
        
        while n < len(ref_seq):
            char_ = ref_seq[n]
            n_gapped = n -gaps
            if (n_gapped >= (offset_start -antisense)):
                if (n_gapped < offset_start):
                    overhang1 += char_
                elif (n_gapped >= (offset_end +sense)):
                    section2 += char_
                elif (n_gapped >= offset_end):
                    overhang2 += char_
                else:
                    highlight_section += char_
            else:
                section1 += char_
        
            if char_ == "-":
                gaps += 1
            n += 1
        
        antisense_overhang_start = len(section1)
        antisense_overhang_end = antisense_overhang_start + len(overhang1)
        color_start_index = antisense_overhang_end
        color_end_index = color_start_index + len(highlight_section)
        sense_overhang_start = color_end_index
        sense_overhang_end = sense_overhang_start + len(overhang2)
        
        print('\na')
        for record in alignment:
            print_highlighted_sequence(record, color_start_index, color_end_index, 
                                       sense_overhang_start, sense_overhang_end, 
                                       antisense_overhang_start, antisense_overhang_end, 
                                       color, overhang_color)

        
def highlight_regions(parser):
    parser.add_option("-b","--bed",action="store", type="string", dest="bed", help="Bed file with genomic coordinates to extract (Required).")
    parser.add_option("-s", "--sense",action="store",type="int",dest="sense",default=0,help="Add a highlighted overhang of this many nucleotides in sense (+) direction of reference strand (Default: 0).")
    parser.add_option("-n", "--antisense",action="store",type="int",dest="antisense",default=0,help="Add a highlighted overhang of this many nucleotides in antisense (-) direction of reference strand (Default: 0).")
    parser.add_option("-c", "--color",action="store",type="string",dest="color",default="GREEN",help="Color for highlighting of annotated regions. Choose from %s (Default: GREEN)." % list(colors_dict.keys()))
    parser.add_option("-g", "--overhang-color",action="store",type="string",dest="overhang", default="",help="Color for highlighting of overhangs. Does nothing if overhangs are length 0 (Default: Same as --color).")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    required = ["bed", ]
    
    for r in required:
        if options.__dict__[r] == None:
            print("You must pass a --%s argument." % r)
            sys.exit()
    
    if options.color not in colors_dict.keys():
        print("Choose color from RED, BLUE, GREEN, YELLOW (Default: GREEN).")
        sys.exit()
    
    if options.overhang == "":
        overhang_color = options.color
    else:
        overhang_color = options.overhang
        if overhang_color not in colors_dict.keys():
            print("Choose overhang color from RED, BLUE, GREEN, YELLOW (Default: Same as --color).")
            sys.exit()
    
    annotation = load_bed_with_range(options.bed)
    alignment_handle = read_maf(handle_)
    
    for alignment in alignment_handle:
        name = alignment[0].id
        df_ = annotation[annotation["sequence"] == name]
        print_highlighted_alignment(alignment, df_, options.sense, options.antisense, options.color, overhang_color)

