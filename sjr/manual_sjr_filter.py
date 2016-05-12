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
good = open('../data/medium_motifs_manual.bed', 'w')

for line in expression:
	chrom, start, end, name, count, strand, expression = line.strip().split('\t')
	start, end, count = int(start), int(end), int(count)

	expression = map(int, expression.split(','))

	e += [(chrom, start, end, strand, expression)]

print 'finsihed reading'


data = open(sys.argv[1], 'r')
for line in data:
	intron = IntronRS(line)
	if not intron.aggt(): continue

	motif = genome_seq[intron.chrom][int(intron.rs) - 30: int(intron.rs) + 30]
	if intron.strand == '-': motif = revcomp(motif)
	motif = motif[10:38]
	motif_score = (score_motif(pwm, motif) - p_min) / (p_max - p_min)
	if not .83 < motif_score <= .87: continue

	coverage = None
	for chrom, start, end, strand, expression in e:
		if chrom  == intron.chrom and start < intron.rs < end and strand == intron.strand:
			coverage = expression
			pos = intron.rs - start if strand == '+' else end - intron.rs
			break

	if coverage != None:
		print motif_score, motif
		if intron.strand == '-':
			coverage = coverage[::-1]
		plt.plot(coverage)
		plt.axvline(pos, c = 'r', linewidth=2)
		plt.show(block = False)

		a = raw_input('any key for yes')

		if a: good.write(str(intron) + '\n')

		plt.close()

	else:
		print motif_score, motif, intron.strand
print c

