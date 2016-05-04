"""
assigns intron by:

1) highest number of splice jxn reads
2) longest intron

Overlooks introns shorted than 1000 bp (what ever is in intron file)
"""

import sys

sjr = open(sys.argv[1], 'r')
introns = open(sys.argv[2], 'r')


anno = {}
for line in introns:
	chrom, start, end, name, juncs, strand = line.strip().split('\t')[:6]
	start, end, juncs = int(start), int(end), int(juncs)
	if strand == '+':
		five, three = start, end
	else:
		five, three = end, start
	key = (chrom, five, strand)

	if key in anno:
		anno[key] += [(three, juncs)]
	else:
		anno[key] = [(three, juncs)]


for line in sjr:
	chrom, start, end, samples, count, strand, seq1, seq2, score1, score2 = line.strip().split('\t')
	start, end, count = int(start), int(end), int(count)
	if strand == '+':
		five, rs = start, end
	else:
		five, rs = end-1, start
	key = (chrom, five, strand)

	if key in anno:
		best = (0, -1)
		for three, juncs in anno[key]:
			if juncs < best[1]: continue
			if strand == '+':
				if rs > three: continue
				if three > best[0] or juncs > best[1]: best = (three, juncs)
			if strand == '-':
				if rs < three: continue
				if three < best[0] or juncs > best[1]: best = (three, juncs)
		if strand == '+':
			begin, stop = five, best[0]
		else:
			begin, stop = best[0], five
		print '\t'.join(map(str, [chrom, begin, stop, samples, rs, strand, count, best[1], seq1, seq2, score1, score2]))
	else:
		# introns less than 1000 not present in express file
		assert abs(five - rs) < 1000
