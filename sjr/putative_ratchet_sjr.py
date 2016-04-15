import pysam
from get_anno_splices import get_jxns
import sys
MIN_INTRON_SIZE = 35

"""
Finds all reads in the given samfile that are broken and have a
5'ss aligning to a known 5'ss and a 3'ss aligning upstream of
a known 3'ss for the given intron.

Constrains that abs(5' - 3') is greater than MIN_INTRON_SIZE.

Need to tune expression for strand based on geometry of library.
"""

samfile = pysam.AlignmentFile(sys.argv[1], 'rb')
sjr = pysam.AlignmentFile(sys.argv[2], 'wb', template = samfile)
bed = open(sys.argv[3], 'w')

ss = get_jxns('../../data/anno.ss')

c = 0
for read in samfile.fetch():
	c += 1
	if not c % 10000000: print c
	blocks = read.get_blocks()
	strand = (read.is_read1 == read.is_reverse)
	if len(blocks) > 1:
		i = 1
		while i < len(blocks) - 1:
			if blocks[i][0] - blocks[i-1][1] < 5:
				blocks[i] = (blocks[i][0], blocks[i+1][1])
				blocks.remove(blocks[i+1])
			else:
				i += 1
		
		for i in xrange(len(blocks) - 1):
			if strand:
				five, three = blocks[i][1] - 1, blocks[i+1][0]
			else:
				three, five = blocks[i][1] - 1, blocks[i+1][0]
			# Some insertion events align to 5'ss and cause false positives

			if abs(five - three) < MIN_INTRON_SIZE: continue


			chrom = 'chr' + samfile.getrname(read.reference_id)
			strand_str = '+' if strand else '-'
			if (chrom, strand_str, five) in ss:
				if three not in ss[(chrom, strand_str, five)]:
					if strand and three < max(ss[(chrom, strand_str, five)]):
						bed.write('\t'.join(map(str, [chrom, five + 1, three, read.query_name, blocks, strand_str])) + '\n')
						sjr.write(read)
					elif not strand and three > min(ss[(chrom, strand_str, five)]):
						bed.write('\t'.join(map(str,[chrom, three + 1, five, read.query_name, blocks, strand_str])) + '\n')
						sjr.write(read)
