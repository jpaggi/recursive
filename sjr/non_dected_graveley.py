import sys

detected = {}

for line in sys.stdin:
	chrom, start, end, sample, entropy, strand, seq1, seq2 = line.strip().split('\t')
	if strand == '+':
		three = int(end) + 1
	else:
		three = int(start)
	detected[(chrom, strand, three)] = entropy



graveley_file = open('../../data/recursiveintrons.bed')
#graveley_file = open('../../data/expressed_graveley_rs_introns.bed', 'r')

for line in graveley_file:
	chrom, start, end, igv, rs, strand = line.strip().split('\t')
	rs = map(int, rs.split(','))
	for i in rs:
		if not (chrom, strand, i) in detected:
			print '\t'.join([chrom, start, end])
