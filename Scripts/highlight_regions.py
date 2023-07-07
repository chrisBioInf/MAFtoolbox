#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:10:29 2023

@author: christopher
"""


import sys
from Bio import AlignIO
import pandas as pd


usage = "\n%prog  [options]"
__version__ = "1.0"

columns = ["sequence", "start", "end", "score", "strand", ]


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def read_maf(filename):
    handle = AlignIO.parse(open(filename, 'r'), format='maf')
    return handle


def load_bed(filename):
    df = pd.read_csv(filename, sep="\t", header=None, names=columns)
    ranges = []
    
    for i in range(0, len(df)):
        start = df["start"].iloc[i]
        end = df["end"].iloc[i]
        ranges.append(set((x for x in range(start, end))))
    
    df["range"] = ranges
    return df 


def print_highlighted_sequence(record, color_start_index, color_end_index):
    seq_ = str(record.seq)
    start = record.annotations.get("start")
    size = record.annotations.get("size")
    strand = record.annotations.get("strand")
    srcsize = record.annotations.get("srcSize")
    section1 = seq_[:color_start_index]
    colored_section = seq_[color_start_index:color_end_index]
    section2 = seq_[color_end_index:]
    print('s %s \t%s \t%s \t%s \t%s \t %s%s%s%s%s' % (record.id, start, size, strand, srcsize,
                                                section1, bcolors.GREEN, 
                                                colored_section, 
                                                bcolors.ENDC, section2))


def print_highlighted_alignment(alignment, df):
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
        colored_section = ''
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
                    colored_section += char_
            else:
                section1 += char_
        
            if char_ == "-":
                gaps += 1
            n += 1
        
        color_start_index = len(section1)
        color_end_index = len(section1) + len(colored_section)
        
        print('\na')
        print_highlighted_sequence(alignment[0], color_start_index, color_end_index)

        for record in alignment[1:]:
            print_highlighted_sequence(record, color_start_index, color_end_index)

        
def highlight_regions(parser):
    parser.add_option("-b","--bed",action="store", type="string", dest="bed", help="Bed file with genomic coordinates to extract (Required).")
    parser.add_option("-m","--maf",action="store", type="string", dest="maf", help="MAF alignment file with coordinates corresponding to bed file (Required).")
    options, args = parser.parse_args()
    
    required = ["bed", "maf"]
    
    for r in required:
        if options.__dict__[r] == None:
            print("You must pass a --%s argument." % r)
            sys.exit()
    
    annotation = load_bed(options.bed)
    handle = read_maf(options.maf)
    
    for alignment in handle:
        name = alignment[0].id
        df_ = annotation[annotation["sequence"] == name]
        print_highlighted_alignment(alignment, df_)
