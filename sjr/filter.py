
from core_pipeline.get_motifs import *
from core_pipeline.load_genome import *

data = open('../data/all_introns.bed', 'r')

genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))

fp_pwm, tp_pwm= make_pwm('../data/anno.ss', genome_seq)

pwm = tp_pwm + fp_pwm

p_min, p_max = get_min_score(pwm), get_max_score(pwm)

for line in data:
	chrom, start, end, samples, rs, strand, sjr, include, exclude, seq1, seq2, score1, score2 = line.strip().split('\t')

	if seq1[-2:] != 'AG' or seq2[:2] != 'GT': continue

	motif = seq1[-20:] + seq2[:8]

	full_score = (score_motif(pwm, motif) - p_min) / (p_max - p_min)

	length = int(end) - int(start)

	samples = samples.split(',')

	recursive_index = float(include) / (float(exclude) + 1)

	print full_score, score1, score2, len(samples), recursive_index, include, length
