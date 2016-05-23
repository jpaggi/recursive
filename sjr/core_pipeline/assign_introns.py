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
		five, three = end+1, start
	key = (chrom, five, strand)

	if key in anno:
		anno[key] += [(three, juncs)]
	else:
		anno[key] = [(three, juncs)]

for line in sjr:
	chrom, start, end, dot, counts, strand = line.strip().split('\t')
	start, end = int(start), int(end)
	if strand == '+':
		five, rs = start, end
	else:
		five, rs = end, start
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
		print '\t'.join(map(str, [chrom, begin, stop, rs, counts, strand, best[1]]))
	else:
		# introns less than 1000 not present in expression file
		assert abs(five - rs) < 1000
