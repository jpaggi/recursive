import pysam
import sys
import matplotlib.pyplot as plt
from math import sqrt, exp

samfile = pysam.AlignmentFile(sys.argv[1])

exons = open(sys.argv[2], 'r')

inserts = []
for line in exons:
	chrom, source, entry, start, end, a, strand = line.strip().split()[:7]
	start, end = int(start), int(end)
	
	for read in samfile.fetch(chrom, start, end):
		if read.is_unmapped: continue
		if read.rname != read.rnext: continue
        if (strand == '+') != (read.is_read1 == read.is_reverse): continue
        if not read.is_proper_pair: continue

        # look for reads straddling jxn
        if strand == '+' and read.is_read2:
        	inner_left = read.get_blocks()[-1][1]
        	inner_right = read.pnext
        	#print '+', inner_right - inner_left
        elif strand == '-' and read.is_read1:
            inner_left = read.pos
            inner_right = read.pnext
            #print '-', inner_right - inner_left
        else:
        	continue

       	#print read.is_read1, read.is_read2, strand

        # print inner_right - inner_left + 100, 'all'

        # print start, inner_left, inner_right, end

        # print abs(read.tlen)
        # print inner_right - inner_left + 100, 'calc'
        #print start, inner_left, inner_right, end
        # if not start < inner_left < end: continue
        # if not start < inner_right < end: continue
        if inner_right - inner_left < 800:	
        	inserts.append(inner_right - inner_left)

        if len(inserts) > 100:
        	break
        print len(inserts)


f = lambda x: 20 * exp(-((x - 126)**2) / (2 * 80.0 ** 2))
print sum(inserts) / float(len(inserts))
print sqrt(sum(map(lambda x: x ** 2, inserts)) / float(len(inserts)) - (sum(inserts) / float(len(inserts))) ** 2)

plt.plot(range(0, 800), map(f, range(0, 800)))

plt.hist(inserts)
plt.show()