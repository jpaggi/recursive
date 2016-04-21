import sys

introns = open(sys.argv[1], 'r')

for line in introns:
	chrom, start, end, name1, name2, strand = line.strip().split('\t')
	start, end = int(start), int(end)

	if end - start > 1000: print line.strip()