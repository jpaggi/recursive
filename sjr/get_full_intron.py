import sys
from get_anno_splices import get_jxns
ss = get_jxns('../../data/anno.ss')


for line in sys.stdin:
	chrom, start, end, sample, offsets, strand, seq1, seq2 = line.strip().split('\t')
	if strand == '+':
		five = int(start) - 1
	else:
		five = int(end)
	print ss[('chr' + chrom, strand, five)]
