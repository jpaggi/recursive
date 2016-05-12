import sys
import matplotlib.pyplot as plt
from core_pipeline.get_motifs import *
from core_pipeline.load_genome import *

def main(rs, plot = True, directory = ''):
	genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
	fp_pwm, tp_pwm= make_pwm('../data/anno.ss', genome_seq)
	pwm = tp_pwm + fp_pwm
	p_min, p_max = get_min_score(pwm), get_max_score(pwm)

	lengths = [0]
	lengths_good_motif = [0]
	lengths_great_motif = [0]
	for intron in rs:
		if not intron.aggt(): continue
		lengths += [intron.length()]
		score = (score_motif(pwm, intron.motif_str(20, 8)) - p_min) / (p_max - p_min)
		if score > .8:
			lengths_good_motif += [intron.length()]
		if score > .9:
			lengths_great_motif += [intron.length()]

	if plot:
		plt.title('Recursively Spliced Intron Lengths')
		plt.xlabel('Intron Length')
		plt.ylabel('Number of Recursive Sites')
		plt.hist(lengths, bins = 50, label = "Total, n = {}".format(len(lengths)))
		plt.hist(lengths_good_motif, color = 'r', bins = 50, label = "Motif > .8, n = {}".format(len(lengths_good_motif)))
		plt.hist(lengths_great_motif, color = 'g', bins = 50, label = "Motif > .9, n = {}".format(len(lengths_great_motif)))
		plt.ylim([0, 130])
		plt.legend()
		if directory:
			plt.savefig(directory+'intron_lengths.png')
		plt.show()

	return lengths, lengths_good_motif, lengths_great_motif
