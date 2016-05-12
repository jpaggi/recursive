import sys

for line in sys.stdin:
	chrom, rs, a, b, c, strand, d, start, end, name, expression, e = line.strip().split()

	if start != '-1':
		print '\t'.join([chrom, start, end, ',', rs, strand, '0', '0', '.', '.', '0', '0'])