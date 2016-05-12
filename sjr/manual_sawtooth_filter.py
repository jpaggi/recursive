import sys
import matplotlib.pyplot as plt
from intron_rs import IntronRS
from core_pipeline.get_motifs import *
from core_pipeline.load_genome import *

grav_file = open('../data/graveley_expression.bed', 'r')
random_scores = map(float, open('../data/full_motif_AGGT.csv', 'r').read().split(',')[:-1])
genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
fp_pwm, tp_pwm= make_pwm('../data/anno.ss', genome_seq)
pwm = tp_pwm + fp_pwm
p_min, p_max = get_min_score(pwm), get_max_score(pwm)

c = 0

expression = open('../data/all_merged.bed', 'r')
e = []
good = open('../data/super_confidence_sawtooth.bed', 'w')

for line in expression:
	chrom, start, end, name, count, strand, expression = line.strip().split('\t')
	start, end, count = int(start), int(end), int(count)

	expression = map(int, expression.split(','))

	e += [(chrom, start, end, strand, expression)]

print 'finsihed reading'


data = open(sys.argv[1], 'r')
for line in data:
	chrom, rs, a, count, score, strand = line.strip().split()
	rs = int(rs)

	if float(score) < .0001: continue
	motif = genome_seq[chrom][int(rs) - 30: int(rs) + 30]
	if strand == '-': motif = revcomp(motif)
	motif = motif[10:38]
	motif_score = (score_motif(pwm, motif) - p_min) / (p_max - p_min)
	if motif_score < .83: continue

	coverage = None
	for e_chrom, e_start, e_end, e_strand, expression in e:
		if chrom  == e_chrom and e_start < rs < e_end and e_strand == strand:
			coverage = expression
			pos = rs - e_start if strand == '+' else e_end - rs
			break

	if coverage != None:
		print motif_score, motif
		if strand == '-':
			coverage = coverage[::-1]
		plt.plot(coverage)
		plt.axvline(pos, c = 'r', linewidth=2)
		plt.show(block = False)

		a = raw_input('any key for yes')

		if a: good.write(line)

		plt.close()

	else:
		print motif_score, motif, intron.strand
print c

