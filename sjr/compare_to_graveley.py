import sys

graveley_file = open('../../data/recursiveintrons.bed')

graveley = []
for line in graveley_file:
	chrom, start, end, igv, rs, strand = line.strip().split('\t')
	rs = map(int, rs.split(','))
	start, end = int(start), int(end)
	graveley += [(chrom, start, end, strand, rs)]


for line in sys.stdin:
	chrom, start, end, sample, entropy, strand, seq1, seq2 = line.strip().split('\t')
	if not int(entropy): continue
	start, end = int(start), int(end)
	overlap = False
	seen = False
	for (g_chrom, g_start, g_end, g_strand, g_rs) in graveley:
		# - strand indexing the same
		# + strand g_rs == end + 1
		if g_chrom == chrom and g_strand == strand and (g_start <= start <= g_end or g_start <= end <= g_end):
			if strand == '+':
				seen = seen or (end + 1 in g_rs)
			else:
				seen = seen or (start in g_rs)
			overlap = True
			# if (end + 1 in g_rs) or (start in g_rs):
			# 	print 'match'
	if not overlap or not seen:
		print chrom, start, end, sample, entropy, seq1, seq2, strand, 'novel'
	else:
		print chrom, start, end, sample, entropy, seq1, seq2, strand, 'graveley'
	# 	pass
	# else:
	# 	print seq1 + seq2
