
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq

from utility import read_maf, check_positional_argument, strand_dict


def write_fasta(seqs, ids, outfile=''):
    records = []
    for (a, b) in zip(seqs, ids):
        if (outfile == ''):
            print('>%s' % b)
            print(a)
        else:
            r = SeqRecord(
                    Seq(a),
                    id=b,
                    name=b,
                    description='',
                    )
            records.append(r)

    if (outfile != ''):
        SeqIO.write(records, open(outfile, 'w'), format='fasta')


def maf_to_fasta(parser):
    parser.add_option("-o","--output",action="store",type="string", default="", dest="outfile",help="FASTA file to write to. If empty, coordinates are redirected to stdout.")
    parser.add_option("-g", "--ungap",action="store_true",default=False,dest="ungap",help="If set, all sequences will be ungapped (Default: False).")
    parser.add_option("-u", "--uppercase",action="store_true",default=False,dest="upper",help="If set, all sequences will be uppercase in output (Default: False).")
    options, args = parser.parse_args() 
    args = args[1:]
    handle_ = check_positional_argument(args)
    alignment_handle = read_maf(handle_)
    seqs = []
    ids = []

    for aln in alignment_handle:
        for record in aln:
            seq = str(record.seq)

            if (options.ungap == True):
                seq = seq.replace('-', '')
            if (options.upper == True):
                seq = seq.upper()

            start = record.annotations.get('start')
            end = start + record.annotations.get('size')
            strand = strand_dict.get(record.annotations.get('strand'))
            id_ = '%s_%s_%s_(%s)' % (
                    str(record.id), 
                    start,
                    end,
                    strand,
                    )
            seqs.append(seq)
            ids.append(id_)

    write_fasta(seqs, ids, options.outfile)

