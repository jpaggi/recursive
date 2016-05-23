import sys
# python sequence/included_exons/merge.py ../data/five_included/motif.bed ../data/ten_included/motif.bed  ../data/twenty_included/motif.bed ../data/total_included/motif.bed
from load_genome import *
from get_motifs import *

"""
adding the monotonicity constraint on the up values doubles the fp rate
but removes a few potentially legitamate values

probably some more sophisticated test using the replicate information that can be used

"""

genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))

fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

class Entry:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom = a[0]
		self.start = int(a[1])
		self.end = int(a[2])
		self.position = a[3]
		self.body = map(int, a[4].split(','))
		self.strand = a[5]
		self.ups = map(int, a[6].split(','))
		self.downs = map(int, a[7].split(','))

for line in sys.stdin:
	entry = Entry(line)

	if sum(entry.ups) < 5: continue

	if entry.position != 'middle': continue

	if sum(entry.downs) * 2 > sum(entry.ups): continue

	diff = [i-j for i, j in zip(entry.body[:-1], entry.body[1:])]
	if not sum(map(lambda x: x < 0, diff)) < 2: continue
	ups_diff = [i-j for i, j in zip(entry.ups[:-1], entry.ups[1:])]
	if not sum(map(lambda x: x < 0, ups_diff)) < 2: continue

	rs = entry.start if entry.strand == '+' else entry.end

	motif = genome_seq[entry.chrom][rs - 20: rs + 20]
	if entry.strand == '-': motif = revcomp(motif)

	motif = motif[:28]

	score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)

	#if score < .8: continue

	print line.strip()
