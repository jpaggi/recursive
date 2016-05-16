import pysam
import sys

samfile = pysam.AlignmentFile(sys.argv[1])

exons = open(sys.argv[2], 'r')

for line in exons:
	chrom, source, entry, start, end, a, strand = line.strip().split()[:7]
	start, end = int(start), int(end)
	for read in samfile.fetch(chrom, start, end):
		if read.is_unmapped: continue
		if read.rname != read.rnext: continue
        if (strand == '+') != (read.is_read1 == read.is_reverse): continue

        # look for reads straddling jxn
        if strand == '+' and read.is_read2:
        	print read.get_blocks()
            inner_left = read.get_blocks()[-1][1]
            inner_right = read.pnext
            #print read.blocks, inner_left, read.pos
        elif strand == '-' and read.is_read1:
            inner_left = read.pos
            inner_right = read.pnext
        else:
        	continue

       	#print read.is_read1, read.is_read2, strand

        # print inner_right - inner_left + 100, 'all'

        # print start, inner_left, inner_right, end

        # print abs(read.tlen)
        # print inner_right - inner_left + 100, 'calc'

        if not start < inner_left < end: continue
        if not start < inner_right < end: continue

        print 'here'


        
        print read.tlen