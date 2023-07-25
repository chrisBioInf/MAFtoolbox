#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 16:23:56 2023

@author: christopher
"""


import pandas as pd

from utility import check_positional_argument


columns = ["match", "mismatch", "rep_match", "Ns", "Q gap count", 
           "Q gap bases", "T gap count", "T gap bases", "strand",
           "Q name", "Q size", "Q start", "Q end", 
           "T name", "T size", "T start", "T end",
           "block count", "blockSizes", "qStarts", "tStarts",
           ]


def load_psl(handle):
    df = pd.read_csv(handle, names=columns, skiprows=5, header=None, sep="\t")
    df["match_fraction"] = df["match"] / df["Q size"]
    return df


def psl_to_bed(parser):
    parser.add_option("-o", "--output", action="store", type="string", dest="output", default="", help="BED file to write to. If empty, results are redirected to stdout.")
    parser.add_option("-m", "--min-match", action="store", type="float", dest="match", default=0, help="Set a minimum fraction of the query sequencce length that must be a match. Other hits will be discarded (Default: 0).")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    seqname = []
    starts = []
    ends = []
    score = []
    strands = []
    count = 0
    
    psl = load_psl(handle_)
    
    for i in range(0, len(psl)):
        if psl["match_fraction"].iloc[i] < options.match:
            continue
        seqname.append(psl["T name"].iloc[i])
        starts.append(psl["T start"].iloc[i])
        ends.append(psl["T end"].iloc[i])
        strands.append(psl["strand"].iloc[i])
        score.append(round(psl["match_fraction"].iloc[i], 2))
        count += 1
    
    if options.output == "":
        for i in range(0, count):
            line = "%s\t%s\t%s\t%s\t%s" % (
                seqname[i], starts[i], ends[i], 
                score[i], strands[i],
                )
            print(line)
    else:
        with open(options.output, 'w') as f:
            for i in range(0, count):
                line = "%s\t%s\t%s\t%s\t%s" % (
                    seqname[i], starts[i], ends[i], 
                    score[i], strands[i],
                    )
                f.write(line + "\n")
