#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 09:03:11 2023

@author: christopher
"""

import pandas as pd
from Bio import AlignIO
import sys


def maf_to_gtf(fs):
    seqname = []
    source = []
    type_ = []
    starts = []
    ends = []
    score = []
    strands = []
    frame = []
    descriptor = []
    filename = fs
    strand_dict = {"forward" : "+", "reverse": "-"}
    
    if "window" in fs:
        this_source = "rnazWindow"
    else:
        this_source = "Multitz"
    
    handle = AlignIO.parse(handle=open(filename, 'r'), format="maf")
    count = 0
    
    for alignment in handle:
        count += 1
        record = alignment[0]
        seqname.append(str(record.id))
        source.append(this_source)
        type_.append("alignment_block")
        starts.append(record.annotations.get("start"))
        ends.append(record.annotations.get("start") + record.annotations.get("size"))
        score.append(".")
        frame.append(".")
        strands.append(record.annotations.get("strand"))
        descriptor.append("Continuously aligned block %s" % count)
        
    
    with open("%s.gtf" % fs, 'w') as f:
        for i in range(0, len(seqname)):
            line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
                seqname[i], source[i], type_[i], starts[i], ends[i], 
                score[i], strands[i], frame[i], descriptor[i]
                )
            print(line)
            f.write(line + "\n")


def main():
    filenames = sys.argv[1:]
    for fs in filenames:
        maf_to_gtf(fs)
    

if __name__ == "__main__":
    main()
