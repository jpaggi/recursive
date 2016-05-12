from matplotlib_venn import venn3
from matplotlib import pyplot as plt
from core_pipeline.get_motifs import *
from core_pipeline.load_genome import *


class Samples:
	def __init__(self, samples = ''):
		self.samples = {}
		for sample in samples:
			kind, time, rep = sample.split('_')
			if time in self.samples:
				if rep in self.samples[time]:
					continue # could add something here later for PEr
				else:
					self.samples[time][rep] = 1
			else:
				self.samples[time] = {rep:1}

	def total_samples(self):
		return sum(self.samples_at_time(time) for time in self.samples)

	def samples_at_time(self, time):
		if time not in self.samples: return 0
		return sum(self.samples[time][rep] for rep in self.samples[time])

	def num_times(self):
		return len(self.samples)

	def union(self, other):
		for time in other.samples:
			if time in self.samples:
				for rep in other.samples[time]:
					if rep in self.samples:
						continue
					else:
						self.samples[time][rep] = 1
			else:
				self.samples[time] = other.samples[time]

def main(data, plot = True, directory = ''):
	grav_file = open('graveley_introns.bed', 'r')
	genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
	fp_pwm, tp_pwm= make_pwm('../data/anno.ss', genome_seq)
	pwm = tp_pwm + fp_pwm
	p_min, p_max = get_min_score(pwm), get_max_score(pwm)

	graveley = {}
	for line in grav_file:
		intron = IntronRS(line)
		graveley[(intron.chrom, intron.strand, intron.rs)] = (1, 0, 0)

	sawtooth_file = open('sawtooth_introns.bed')
	sawtooth = {}
	for line in sawtooth_file:
		intron = IntronRS(line)
		key = (intron.chrom, intron.strand, intron.rs)

		if key in graveley:
			graveley[key] = (1, 1, 0)
			sawtooth[key] = (1, 1, 0)
		else:
			sawtooth[key] = (0, 1, 0)

	detected = {}
	for intron in data:
		key = (intron.chrom, intron.strand, intron.rs)
		if key in graveley:
			a, b, c = graveley[(intron.chrom, intron.strand, intron.rs)]
			graveley[key] = (a, b, 1)
			detected[key] = (a, b, 1)
		if key in sawtooth:
			a, b, c = sawtooth[(intron.chrom, intron.strand, intron.rs)]
			sawtooth[key] = (a, b, 1)
			detected[key] = (a, b, 1)
		if key not in graveley and key not in sawtooth:
			detected[key] = (0, 0, 1)

	All = len(filter(lambda x: x == (1, 1, 1), graveley.values()))
	print 'All            ', All
	GS = len(filter(lambda x: x == (1, 1, 0), graveley.values()))
	print 'Grav + sawtooth', GS
	GSJR = len(filter(lambda x: x == (1, 0, 1), graveley.values()))
	print 'Grav + SJR     ', GSJR
	SSJR = len(filter(lambda x: x == (0, 1, 1), sawtooth.values()))
	print 'SJR + sawtooth ', SSJR

	G = len(filter(lambda x: x == (1, 0, 0), graveley.values()))
	print 'Grav Only      ', G
	SJR = len(filter(lambda x: x == (0, 0, 1), detected.values()))
	print 'SJR Only       ', SJR
	S = len(filter(lambda x: x == (0, 1, 0), sawtooth.values()))
	print 'Sawtooth Only  ', S

	print filter(lambda x: sawtooth[x] == (1, 1, 0), sawtooth.keys())

	G += -1

	venn3(subsets = (G, SJR, GSJR, S, GS, SSJR, All), set_labels = ('Graveley', 'SJR', 'Sawtooth'))

	plt.savefig('../data/venn/manually_filtered_high_expression.png')

	plt.show()

if __name__ == "__main__":
	import sys
	from intron_rs import IntronRS
	data = open(sys.argv[1], 'r')
	rs = [IntronRS(line) for line in data]
	print main(rs, True)
