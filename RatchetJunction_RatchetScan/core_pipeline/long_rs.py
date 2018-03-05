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

		#if abs(int_start - start) < 5: print int_start - start, strand
		if strand == '+' and (int_start + 1000 > start or int_end != end-1): continue
		if strand == '-' and (int_start != start+1 or end > int_end - 1000): continue
		print line.strip()
		break