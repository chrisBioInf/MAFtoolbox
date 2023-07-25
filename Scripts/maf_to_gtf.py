#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 09:03:11 2023

@author: christopher
"""


from maf_to_bed import coordinate_overlap
from utility import read_maf, check_positional_argument, strand_dict


def maf_to_gtf(fs, output, cluster):
    seqname = []
    source = []
    type_ = []
    starts = []
    ends = []
    score = []
    strands = []
    frame = []
    descriptor = []

    this_source = "MAF"
    
    alignment_handle = read_maf(fs)
    count = 0
    
    if cluster == False:
        for alignment in alignment_handle:
            count += 1
            record = alignment[0]
            seqname.append(record.id)
            source.append(this_source)
            type_.append("alignment_block")
            starts.append(record.annotations.get("start"))
            ends.append(record.annotations.get("start") + record.annotations.get("size"))
            score.append(".")
            frame.append(".")
            strands.append(strand_dict.get(record.annotations.get("strand")))
            descriptor.append("Continuously aligned block %s" % count)
            
    else:
        alignment = next(alignment_handle)
        record = alignment[0]
        seqname.append(record.id)
        type_.append("alignment_block")
        source.append(this_source)
        starts.append(record.annotations.get("start"))
        ends.append(record.annotations.get("start") + record.annotations.get("size"))
        score.append(".")
        frame.append(".")
        strands.append(strand_dict.get(record.annotations.get("strand")))
        descriptor.append("Continuously aligned block %s" % count)
        count += 1
        
        for alignment in alignment_handle:
            record = alignment[0]
            start = record.annotations.get("start")
            end = record.annotations.get("start") + record.annotations.get("size")
            strand = strand_dict.get(record.annotations.get("strand"))
            
            if (record.id == seqname[-1]) and (strand == strands[-1]) and coordinate_overlap(ends[-1], start, end):
                ends[-1] = end
                continue
            
            seqname.append(record.id)
            type_.append("alignment_block")
            source.append(this_source)
            starts.append(start)
            ends.append(end)
            frame.append(".")
            score.append(".")
            strands.append(strand)
            descriptor.append("Continuously aligned block %s" % count)
            count += 1
        
    if output == "":
        for i in range(0, count):
            line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
                seqname[i], source[i], type_[i], starts[i], ends[i], 
                score[i], strands[i], frame[i], descriptor[i]
                )
            print(line)
    else:
        with open(output, 'w') as f:
            for i in range(0, count):
                line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
                    seqname[i], source[i], type_[i], starts[i], ends[i], 
                    score[i], strands[i], frame[i], descriptor[i]
                    )
                f.write(line + "\n")


def write_to_gtf(parser):
    parser.add_option("-o","--output",action="store",type="string", default="", dest="out_file",help="GTF file to write to. If empty, coordinates are redirected to stdout.")
    parser.add_option("-c", "--cluster",action="store_true",default=False,dest="cluster",help="If set, overlapping and same-stranded coordinates will be merged (Default: False).")
    options, args = parser.parse_args()
    args = args[1:]
    handle_ = check_positional_argument(args)
    
    maf_to_gtf(handle_, output=options.out_file, cluster=options.cluster)

