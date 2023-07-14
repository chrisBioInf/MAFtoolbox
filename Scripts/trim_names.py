#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:13:24 2023

@author: christopher
"""


from copy import deepcopy

from utility import read_maf, write_maf, print_maf_alignment, check_positional_argument


def trim_ids(parser):
    parser.add_option("-o","--output",action="store",type="string", default="", dest="out_file",help="MAF file to write to. If empty, results alignments are redirected to stdout.")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    alignment_handle = read_maf(handle_)
    records = []
    
    for alignment in alignment_handle:
        record = deepcopy(alignment)
        name = str(alignment[0].id).split(".")[1] + ".1"
        record[0].id = name
        
        if options.out_file == "":
            print_maf_alignment(record)
        else:
            records.append(record)
    
    if len(records)> 0:
        write_maf(records, options.out_file,)
    

