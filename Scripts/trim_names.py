#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:13:24 2023

@author: christopher
"""


from copy import deepcopy
from Bio import AlignIO
import sys


def trim_this(filename):
    handle = AlignIO.parse(open(filename, 'r'), format='maf')
    records = []
    
    for alignment in handle:
        record = deepcopy(alignment)
        record[0].id = str(alignment[0].id).split(".")[1] + ".1"
        records.append(record)

    AlignIO.write(records, open(filename, 'w'), format='maf')
    

trim_this(sys.argv[1])
