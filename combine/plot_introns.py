from intron import Intron
import matplotlib.pyplot as plt
from load_genome import *
from get_motifs import *


def main(entries):
	genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
	fp_pwm, tp_pwm = make_pwm('../data/annotations/anno.ss', genome_seq)
	pwm = tp_pwm + fp_pwm
	min_score, max_score = get_min_score(pwm), get_max_score(pwm)


	introns = []
	for line in open('../data/coverage/expressed_coverage.bed'):
		chrom, start, end, name, counts, strand, coverage = line.strip().split('\t')
		coverage = map(int, coverage.split(','))
		introns += [Intron(chrom, int(start), int(end), strand, int(counts), coverage)]

	introns = sorted(introns, key = lambda x: x.length(), reverse = True)

	for entry in entries:
		chrom, rs, strand = entry.chrom, entry.rs, entry.strand
		for intron in introns:
			if intron.compatible(chrom, strand, rs):
				exon = sum(entry.down_counts) * 5 > sum(entry.junc) 
				intron.add(rs, exon)



	for intron in introns:
		if intron.chrom != '3L' or intron.start != 12021978: continue
		print intron.chrom, intron.start, intron.end
		seq = genome_seq[intron.chrom][intron.start:intron.end]
		if intron.strand == '-':
			intron.coverage = intron.coverage[::-1]
			seq = revcomp(seq)

		plt.plot(intron.coverage)

		for rs, exon in zip(intron.rs, intron.exon):
			if intron.strand == '+':
				pos = rs - intron.start
			else:
				pos = intron.end - rs
			c = 'r' if not exon else 'g'
			plt.axvline(pos, linewidth = 2, c = c)

		for i in range(30, len(seq) - 30):
			motif = seq[i - len(tp_pwm): i + len(fp_pwm)]
			if motif[len(tp_pwm)-2: len(tp_pwm)+2] != 'AGGT': continue
			score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)
			if score > .8: plt.scatter([i], [1], marker = 'o', linewidths = [(score - .8) * 100], c = 'k')
		plt.autoscale(tight = True)
		plt.ylim([0, 20])
		plt.show(block = False)
		a = raw_input("enter to continue")

		plt.savefig('../data/3L:12021978.png')
		plt.close()


if __name__ == '__main__':
	import sys
	from standard_table_reader import Entry
	main([Entry(line) for line in open(sys.argv[1], 'r')])
