#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:31:30 2023

@author: christopher
"""


from Bio import AlignIO
import sys

from utility import load_bed, write_maf, print_maf_alignment


def extract_blocks(maf_file, annotations, output):
    extracted_blocks = []
    handle = AlignIO.parse(open(maf_file, 'r'), format="maf")
    annotations["annotation_key"] = annotations['sequence'] + '_' + annotations['start'].astype(str) + '_' + annotations['end'].astype(str) 
    
    for alignment in handle:
        ref_seq = alignment[0]
        name = str(ref_seq.id)
        start = int(ref_seq.annotations["start"])
        end = start + int(ref_seq.annotations["size"])
        key = "%s_%s_%s" % (name, start, end)
        
        if key in list(annotations['annotation_key']):
            if output == "":
                print_maf_alignment(alignment)
            else:
                extracted_blocks.append(alignment)
            
    return extracted_blocks


def extract_alignment(parser):
    if len(sys.argv) < 2:
        print("Usage:\npython extract_maf_by_coordinate.py -b [BED FILE] -a [ALIGNMENT FILE] ")
        sys.exit()
    parser.add_option("-b","--bed",action="store", type="string", dest="bed", help="Bed file with genomic coordinates to extract (Required).")
    parser.add_option("-m","--maf",action="store", type="string", dest="maf", help="MAF alignment file with coordinates corresponding to bed file (Required).")
    parser.add_option("-o","--output",action="store",type="string", default="", dest="output", help="MAF file to write to. If empty, results alignments are redirected to stdout.")
    options, args = parser.parse_args()
    
    required = ["bed", "maf"]
    
    for r in required:
        if options.__dict__[r] == None:
            print("You must pass a --%s argument." % r)
            sys.exit()
            
    annotations = load_bed(options.bed)
    extracted_blocks = extract_blocks(options.maf, annotations, options.output)
    
    if len(extracted_blocks) > 0:
        write_maf(extracted_blocks, options.output)
    
