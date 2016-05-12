import sys


graveley = open('../data/graveley.bed', 'r')

for line in graveley:
	chrom, start, end, name, count, strand = line.strip().split()
	rs = int(end) if strand == '+' else int(start)

	print '\t'.join(map(str, [chrom, rs, rs + 1, name, count, strand]))
