#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:31:30 2023

@author: christopher
"""


from optparse import OptionParser
import pandas as pd
from Bio import AlignIO
import sys


__version__ = "0.5"

columns = ["sequence", "start", "end", "score", "strand", ]
columns_2 = ["sequence", "source", "type", "start", "end", "score", "empty", "strand", "attribute"]


def load_bed(filename):
    try:
        df = pd.read_csv(filename, sep="\t", header=None, names=columns_2)
    except Exception:
        df = pd.read_csv(filename, sep="\t", header=None, names=columns)
    return df


def extract_blocks(maf_file, annotations):
    
    extracted_blocks = []
    handle = AlignIO.parse(open(maf_file, 'r'), format="maf")
    annotation_index = 0
    
    annotation_name = annotations["sequence"].iloc[annotation_index]
    annotation_start = annotations["start"].iloc[annotation_index]
    annotation_end = annotations["end"].iloc[annotation_index]
    
    for alignment in handle:
        ref_seq = alignment[0]
        name = str(ref_seq.id)
        # alignment[0].id = name
        start = int(ref_seq.annotations["start"])
        end = start + int(ref_seq.annotations["size"])
        
        if (name == annotation_name) and (start == annotation_start) and (end == annotation_end):
            print(alignment)
            extracted_blocks.append(alignment)
            annotation_index += 1
            annotation_name = annotations["sequence"].iloc[annotation_index]
            annotation_start = annotations["start"].iloc[annotation_index]
            annotation_end = annotations["end"].iloc[annotation_index]
            
            if annotation_index > len(annotations):
                break
    
    return extracted_blocks


def main():
    if len(sys.argv) < 2:
        print("Usage:\npython extract_maf_by_coordinate.py -b [BED FILE] -a [ALIGNMENT FILE] ")
        sys.exit()
    usage = "\npython %prog  -b [BED FILE] -a [ALIGNMENT FILE] \nNote that program assumes both files to be ordered."
    parser = OptionParser(usage, version="%prog " + __version__)
    parser.add_option("-b","--bed",action="store", type="string", dest="bed", help="Bed file with genomic coordinates to extract.")
    parser.add_option("-a","--maf",action="store", type="string", dest="maf", help="MAF alignment file with coordinates corresponding to bed file.")
    parser.add_option("-o","--output",action="store",type="string", dest="output", help="Designate an output file.")
    
    required = ["bed", "maf", "output"]
    options, args = parser.parse_args()
    
    for r in required:
        if options.__dict__[r] == None:
            print("You must pass a --%s argument." % r)
            sys.exit()
            
    annotations = load_bed(options.bed)
    extracted_blocks = extract_blocks(options.maf, annotations)
    
    AlignIO.write(extracted_blocks, open(options.output, "w"), format="maf")


if __name__ == "__main__":
    main()
