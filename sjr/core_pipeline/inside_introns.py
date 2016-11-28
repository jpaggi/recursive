import sys

introns = {'+': {}, '-': {}}

for line in open(sys.argv[1]):
	chrom, start, end, a, b, strand = line.strip().split()[:6]
	start, end = int(start), int(end)
	if not chrom in introns[strand]: introns[strand][chrom] = []
	introns[strand][chrom] += [(start, end)]

for line in sys.stdin:
	chrom, start, end, a, b, strand = line.strip().split('\t')[:6]
	start, end = int(start), int(end)
	if not chrom in introns[strand]: continue
	for int_start, int_end in introns[strand][chrom]:
		if int_start + 2 < start < int_end - 2 and int_start + 2 < end < int_end - 2:
			print line.strip()