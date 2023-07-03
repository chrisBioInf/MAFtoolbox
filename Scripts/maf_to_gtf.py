#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 09:03:11 2023

@author: christopher
"""

import pandas as pd
from optparse import OptionParser
from Bio import AlignIO
import sys



usage = "\n%prog  [options]"
__version__ = "1.0"

def maf_to_gtf(fs, trim_name=False):
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
        this_source = "MAF"
    
    handle = AlignIO.parse(handle=open(filename, 'r'), format="maf")
    count = 0
    
    for alignment in handle:
        count += 1
        record = alignment[0]
        
        if trim_name == True:
            seqname.append(str(record.id).split(".")[1] + ".1")
        else:
            seqname.append(str(record.id))
        source.append(this_source)
        type_.append("alignment_block")
        starts.append(record.annotations.get("start"))
        ends.append(record.annotations.get("start") + record.annotations.get("size"))
        score.append(".")
        frame.append(".")
        strands.append(record.annotations.get("strand"))
        descriptor.append("Continuously aligned block %s" % count)
        
    
    if fs.endswith(".maf"):
        outfile = fs.replace(".maf", ".gtf")
    else:
        outfile = fs + ".gtf"
    
    with open(outfile, 'w') as f:
        for i in range(0, len(seqname)):
            line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
                seqname[i], source[i], type_[i], starts[i], ends[i], 
                score[i], strands[i], frame[i], descriptor[i]
                )
            f.write(line + "\n")


def main():
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--input",action="store",type="string", dest="in_file",help="The (MAF) input file (Required).")
    parser.add_option("-t", "--trim", action="store", default=False, dest="trim", help="If sequence names should be split by seperator. Default: false")
    options, args = parser.parse_args()
    
    maf_to_gtf(options.in_file, trim_name=options.trim)
    

if __name__ == "__main__":
    main()
