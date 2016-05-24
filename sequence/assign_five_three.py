"""
assigns intron by:

1) highest number of splice jxn reads
2) longest intron

Overlooks introns shorted than 1000 bp (what ever is in intron file)
"""

import sys
from standard_table_reader import Entry

entries = open(sys.argv[1], 'r')
introns = open('../data/all_merged.bed', 'r')

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

for line in entries:
	entry = Entry(line)
	chrom, five, three, rs, strand = entry.chrom, entry.five(), entry.three(), entry.rs, entry.strand

	if five == set([0]):
		#assign five and three
		best = (-1, -1, -1)
		for key in anno:
			a_chrom, a_five, a_strand = key
			if a_chrom != chrom or a_strand != strand: continue
			for a_three, a_juncs in anno[key]:
				if a_juncs < best[2]: continue
				if strand == '+' and a_five + 5 < rs < a_three - 5:
					best = (a_five, a_three, a_juncs)
				elif strand == '-' and a_three + 5 < rs < a_five - 5:
					best = (a_five, a_three, a_juncs)

		if best[2] != -1:
			entry.set_five(best[0])
			entry.set_three(best[1])
			entry.set_intron_sjr(best[2])
			#print entry
	else:
		# assign three
		best = (0, -1)
		for f in five:
			key = (chrom, f, strand)
			if key in anno:
				for three, juncs in anno[key]:
					if juncs < best[1]: continue
					if strand == '+':
						if rs >= three + 5: continue
						if three > best[0] or juncs > best[1]: best = (three, juncs)
					if strand == '-':
						if rs + 5 <= three: continue
						if three < best[0] or juncs > best[1]: best = (three, juncs)
		if best[1] != -1:
			entry.set_three(best[0])
			#print entry
		else:	
			# introns less than 1000 not present in expression file
			print entry

			print [(chrom, f + i, strand) in anno for i in range(-2, 3)]
			print f, rs
			#assert abs(f - rs) < 1000
