#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:09:33 2023

@author: christopher
"""


import sys
import os
from optparse import OptionParser

__author__ = "Christopher Klapproth"
__institution__= "University Leipzig"
__credits__ = []
__license__ = "GPLv2"
__version__="0.5.0"
__maintainer__ = "Christopher Klapproth"
__email__ = "christopher@bioinf.uni-leipzig.de"
__status__ = "Development"


THIS_FILE = os.path.abspath(__file__)
THIS_DIR = os.path.dirname(THIS_FILE)
SCRIPT_DIR = os.path.join(THIS_DIR, "Scripts")

if SCRIPT_DIR not in sys.path:
    sys.path.append(SCRIPT_DIR)
    

from wga_block_merger import merge
from removeRepeats import mask_repeat_regions
from select_maf_sequences import select_seqs
from maf_to_gtf import write_to_gtf
from maf_to_bed import write_to_bed
from extract_maf_by_coordinates import extract_alignment
from highlight_regions import highlight_regions
from trim_names import trim_ids    


usage_statement = "Usage: MAFtools [program] [options], with program being one of 'blockmerge', 'mask', 'filter', 'extract', 'highlight', 'toBed', 'toGTF', 'trimNames'."


def main():
    parser = OptionParser(usage=usage_statement, version="__version__")
    args = sys.argv
    
    if len(args) < 2:
        print(usage_statement)
        sys.exit()
    if args[1] == "blockmerge":
        merge(parser)
    elif args[1] == "mask":
        mask_repeat_regions(parser)
    elif args[1] == "filter":
        select_seqs(parser)
    elif args[1] == "extract":
        extract_alignment(parser)
    elif args[1] == "highlight":
        highlight_regions(parser)
    elif args[1] == "toBed":
        write_to_bed(parser)
    elif args[1] == "toGTF":
        write_to_gtf(parser)
    elif args[1] == "trimNames":
        trim_ids(parser)
    else:
        print(usage_statement)
        sys.exit()
    

if __name__=='__main__':
    main()
