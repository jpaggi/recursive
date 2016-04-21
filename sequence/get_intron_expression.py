import sys
import pysam

"""
Currently asserting that all introns are much 
longer than the insert length
"""

MAX_INSERT = 300
INSIDE = True
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
introns =  open(sys.argv[2], 'r')

def merge_blocks(blocks):
    """
    Removes short gaps in read that are likely indels
    """
    i = 1
    while i < len(blocks) - 1:
        if blocks[i][0] - blocks[i-1][1] < 5:
            blocks[i] = (blocks[i][0], blocks[i+1][1])
            blocks.remove(blocks[i+1])
        else:
            i += 1
    return blocks

for intron in introns:
    chrom, intron_start, intron_end, name1, name2, strand = intron.strip().split('\t')
    intron_start = int(intron_start)
    intron_end = int(intron_end)
    expression = [0] * (intron_end - intron_start)
    sjr_count = 0
    for read in samfile.fetch(chrom, intron_start-MAX_INSERT, intron_end+MAX_INSERT):
        # look for reads straddling jxn
        if strand == '+' and read.is_read2:
            if read.rname != read.rnext: continue
            inner_left = read.blocks[-1][1]
            inner_right = read.pnext

            consistent = inner_left + 5 < intron_start and intron_end < inner_right - 5
            if consistent and intron_start - inner_left + inner_right - intron_end < MAX_INSERT:
                assert inner_right - inner_left > 2 * MAX_INSERT
                sjr_count += 1
        elif strand == '-' and read.is_read1:
            if read.rname != read.rnext: continue
            inner_left = read.pos
            inner_right = read.pnext

            consistent = inner_left + 5 < intron_start and intron_end < inner_right - 5
            if consistent and intron_start - inner_left + inner_right - intron_end < MAX_INSERT:
                assert inner_right - inner_left > 2 * MAX_INSERT
                sjr_count += 1

        # look for split reads spanning intron
        blocks = merge_blocks(read.get_blocks())
        for i in range(len(blocks) - 1):
            if blocks[i][1] == intron_start and blocks[i+1][0] == intron_end + 1:
                sjr_count += 1

        # add expression between two read ends
        if INSIDE and read.is_proper_pair and read.is_read2:
            if strand == '+':
                inner_left = read.blocks[-1][1]
                inner_right = read.pnext
            else:
                inner_left  = read.pnext
                inner_right = read.pos
            start = max(intron_start, inner_left)
            end = min(intron_end, inner_right)
            if inner_right - inner_left < MAX_INSERT:
                for i in xrange(start, end):
                    expression[i - intron_start] += 1

        # add expression for current read end
        for block in blocks:
            start = max(intron_start, block[0])
            end = min(intron_end, block[1])
            for i in xrange(start, end):
                expression[i - intron_start] += 1

    print '\t'.join(map(str, [chrom, intron_start, intron_end, name1, sjr_count, strand]) + [','.join(map(str, expression))])
