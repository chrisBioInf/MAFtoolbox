#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 11:32:08 2023

@author: christopher
"""

import pandas as pd
import sys


names = ["chromosome", "source", "type", "start", "end", "score", "strand", "empty", "locus"]


def merge_overlapping_blocks(df, filename):
    chrs = []
    starts = []
    ends = []
    ids = []
    scores = []
    strands = []
    
    def in_range(end, start1, end1):
        if (start1 < end) and (end < end1):
            return True
        return False
    
    def find_longest_window(end, df_, offset):
        n = 0
        new_end = end
        while n < len(df_):
            if in_range(end, df_["start"].iloc[n], df_["end"].iloc[n]):
                new_end = df_["end"].iloc[n]
                n += 1
            else:
                return new_end, n
        return new_end, n
    
    i, j = 0, 0
    offset = 1
    count = 0
    
    while (i < len(df)-1) and (j < len(df)-1):
        j = i+1
        score = df["score"].iloc[i]
        end = df["end"].iloc[i]
        new_end, n = find_longest_window(end, df[j:], offset)
        count += 1
        
        chrs.append(df["chromosome"].iloc[i])
        starts.append(df["start"].iloc[i])
        ends.append(new_end)
        ids.append("locus%s" % count)
        strands.append(".")
        scores.append(score)
        i = i+n+1
            
    filename = filename.replace(".gtf", ".bed")
    
    with open("%s" % filename, 'w') as f:
        for i in range(0, len(chrs)):
            line = "%s\t%s\t%s\t%s\t%s\t%s" % (
                chrs[i], starts[i], ends[i], ids[i], scores[i], strands[i],
                )
            f.write(line + "\n")


def main():
    filenames = sys.argv[1:]
    for fs in filenames:
        df = pd.read_csv(fs, sep="\t", names=names, index_col=None)
        # df.sort_values(by="start", inplace=True, ascending=True)
        merge_overlapping_blocks(df, fs)
    

if __name__ == "__main__":
    main()
