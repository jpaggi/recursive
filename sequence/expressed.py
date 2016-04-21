import sys

for line in sys.stdin:
	chrom, start, end, name, sjr, strand, expression = line.strip().split('\t')

	if int(sjr) and int(end) - int(start) > 1000:
		print line.strip()
