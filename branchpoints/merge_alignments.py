from sam_iter import SamIter
import sys

def unique_aligns(aligns):
    out = []
    for align in aligns:
        if align not in out:
            out.append(align)
    return out

def format_entry(qname, seq, qual, aligns):
    aligns = unique_aligns(aligns)
    return '\t'.join([cur_qname, cur_seq, cur_qual, ' '.join(aligns)])

def complement(char):
    if char == 'A': return 'T'
    if char == 'C': return 'G'
    if char == 'G': return 'C'
    if char == 'T': return 'A'
    return char

def revcomp(seq):
    return ''.join(map(complement, seq[::-1]))

cur_qname = ''
aligns = []
cur_seq = ''
cur_qual = ''
for read in SamIter(sys.argv[1], unmapped = True):
    if read.qname == cur_qname:
        aligns += [read.align_str()]
    else:
        if cur_qname:
            print format_entry(cur_qname, cur_seq, cur_qual, aligns)
            
        cur_qname = read.qname
        cur_seq = revcomp(read.seq) if read.is_reverse else read.seq
        cur_qual = read.qual[::-1] if read.is_reverse else read.qual
        aligns = [read.align_str()]

print format_entry(cur_qname, cur_seq, cur_qual, aligns)
