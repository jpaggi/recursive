import sys
import matplotlib.pyplot as plt
from intron_rs import IntronRS
from core_pipeline.get_motifs import *
from core_pipeline.load_genome import *


def get_cdf(ri):
	x_vals = [.01 * x for x in range(101)]
	y_vals = []
	for x in x_vals:
		y_vals += [len(filter(lambda z: z[0] > x, ri))]
	return x_vals, y_vals

def main(rs, plot = True, directory = ''):

	genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
	fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)
	pwm = tp_pwm + fp_pwm
	min_score, max_score = get_min_score(pwm), get_max_score(pwm)

	grav_ris = []
	graveley = {}
	grav_file = open('../data/graveley.bed', 'r')
	for line in grav_file:
		chrom, start, end, name, score, strand = line.strip().split('\t')
		ratchet = int(start) if strand == '-' else int(end)
		graveley[(chrom, strand, ratchet)] = (set(), 0, '', '')

	ri = []
	for intron in rs:
		if intron.expression() == 0: continue
		if not intron.aggt(): continue
		ri += [(intron.recursive_index(), (score_motif(pwm, intron.motif_str(20, 8)) - min_score) / (max_score - min_score))]
		# if intron.recursive_index() > .1:
		# 	print intron.motif_str()

		key = (intron.chrom, intron.strand, intron.rs)
		if key in graveley:
			grav_ris += [ri[-1]]


	plt.title('Recursive Index and Motif Strength')
	plt.xlabel('Recursive Index')
	plt.ylabel('Number of Recursive Sites')

	for motif_score, color in zip([0, .7, .8, .85, .9], ['r', 'y', 'c', 'b', 'm']):
		above_thresh = filter(lambda z: z[1] >= motif_score, ri)
		x_vals, y_vals = get_cdf(above_thresh)

		plt.plot(x_vals, y_vals, color = color, label = "{}, n = {}".format(motif_score, len(above_thresh)))

	x_vals, y_vals = get_cdf(grav_ris)
	plt.plot(x_vals, y_vals, color = 'g', label = "{}, n = {}".format('Graveley (Any Motif)', len(grav_ris)))

	plt.ylim([0, 300])
	plt.legend(title = 'Motif Score Cutoff')

	if directory:
		plt.savefig(directory+'recursive_index.png')
	plt.show()
