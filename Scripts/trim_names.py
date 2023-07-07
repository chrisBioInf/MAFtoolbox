#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:13:24 2023

@author: christopher
"""


import sys
from copy import deepcopy

from utility import read_maf, write_maf, print_maf_alignment


def trim_ids(parser):
    parser.add_option("-i","--input",action="store",type="string", dest="input",help="The (MAF) input file (Required).")
    parser.add_option("-o","--output",action="store",type="string", default="", dest="out_file",help="MAF file to write to. If empty, results alignments are redirected to stdout.")
    options, args = parser.parse_args()
    
    required = ["input"]
    
    for r in required:
        if options.__dict__[r] == None:
            print("You must pass a --%s argument." % r)
            sys.exit()
    
    handle = read_maf(options.input)
    records = []
    
    for alignment in handle:
        record = deepcopy(alignment)
        name = str(alignment[0].id).split(".")[1] + ".1"
        record[0].id = name
        
        if options.out_file == "":
            print_maf_alignment(record)
        else:
            records.append(record)
    
    if len(records)> 0:
        write_maf(records, options.out_file,)
    

