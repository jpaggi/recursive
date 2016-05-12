import matplotlib.pyplot as plt
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

def plot_set(counts, times, label, timepoints, color):
	y_vals = []
	x_vals = []
	ticks = []
	for i, time in enumerate(timepoints):
		begin = i * 4 + 1
		y_vals += counts[time]
		ticks += map(str, range(1, len(counts[time]) + 1))
		x_vals += range(begin, begin + len(counts[time]))

	begin = (len(timepoints) + 1) * 3 + 1
	x_vals += range(begin, begin + len(times))
	y_vals += times
	ticks += map(str, range(5))
	plt.bar(x_vals, y_vals, align = 'center', color = color, tick_label = ticks, label = label)

def main(data, plot = True, directory = ''):
	grav_file = open('../data/graveley_expression.bed', 'r')
	random_scores = map(float, open('../data/full_motif_AGGT.csv', 'r').read().split(',')[:-1])
	genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
	fp_pwm, tp_pwm= make_pwm('../data/anno.ss', genome_seq)
	pwm = tp_pwm + fp_pwm
	p_min, p_max = get_min_score(pwm), get_max_score(pwm)

	graveley = {}
	for line in grav_file:
		chrom, start, end, name, score, strand = line.strip().split('\t')
		rs = int(start) if strand == '-' else int(end)
		graveley[(chrom, strand, rs)] = (Samples(), 0, '', '', score)


	detected = {}
	for intron in data:
		if not intron.aggt(): continue
		if (intron.chrom, intron.strand, intron.rs) in graveley:
			samples, total, a, b, score = graveley[(intron.chrom, intron.strand, intron.rs)]
			samples.union(Samples(intron.samples))
			graveley[(intron.chrom, intron.strand, intron.rs)] = (samples, total + intron.sjr, intron.seq1, intron.seq2, score)

		if (intron.chrom, intron.strand, intron.rs) in detected:
			sample, total, seq1, seq2 = detected[(intron.chrom, intron.strand, intron.rs)]
			samples.union(Samples(intron.samples))
			detected[(intron.chrom, intron.strand, intron.rs)] = (samples, total + intron.sjr, intron.seq1, intron.seq2)
		else:
			samples = Samples(intron.samples)
			detected[(intron.chrom, intron.strand, intron.rs)] = (samples, intron.sjr, intron.seq1, intron.seq2)


	timepoints = ['5', '10', '20', 'total']

	total_counts = {}
	total_times = [0] * (len(timepoints) + 1)
	novel_scores_times = [[],[],[],[]]
	novel_scores = {}
	for time in timepoints:
		total_counts[time] = [0] * 3
		novel_scores[time] = [[],[],[]]
	total_counts['total'] = [0] * 2
	novel_scores['total'] = [[],[]]

	for key in detected:
		samples, total, seq1, seq2 = detected[key]
		score = (score_motif(pwm, seq1[-20:] + seq2[:8]) - p_min) / (p_max - p_min)
		for time in timepoints:
			count = samples.samples_at_time(time)
			if count:
				total_counts[time][count-1] += 1
				if key not in graveley: novel_scores[time][count-1] += [score]
		total_times[samples.num_times()] += 1
		if key not in graveley: novel_scores_times[samples.num_times()-1] += [score]



	motif_counts = {}
	motif_times = [0] * (len(timepoints) + 1)
	for time in timepoints:
		motif_counts[time] = [0] * 3
	motif_counts['total'] = [0] * 2

	for key in detected:
		samples, total, seq1, seq2 = detected[key]
		score = (score_motif(pwm, seq1[-20:] + seq2[:8]) - p_min) / (p_max - p_min)
		if score < .87: continue
		for time in timepoints:
			count = samples.samples_at_time(time)
			if count: 
				motif_counts[time][count-1] += 1			
		motif_times[samples.num_times()] += 1
		
	weak_motif_counts = {}
	weak_motif_times = [0] * (len(timepoints) + 1)
	for time in timepoints:
		weak_motif_counts[time] = [0] * 3
	weak_motif_counts['total'] = [0] * 2

	for key in detected:
		samples, total, seq1, seq2 = detected[key]
		score = (score_motif(pwm, seq1[-20:] + seq2[:8]) - p_min) / (p_max - p_min)
		if score < .8: continue
		for time in timepoints:
			count = samples.samples_at_time(time)
			if count:
				weak_motif_counts[time][count-1] += 1
		weak_motif_times[samples.num_times()] += 1

	grav_counts = {}
	grav_times = [0] * (len(timepoints) + 1)
	grav_scores = []
	for time in timepoints:
		grav_counts[time] = [0] * 3
	grav_counts['total'] = [0] * 2

	for key in graveley:
		samples, total, seq1, seq2, sjr_count = graveley[key]
		#if int(sjr_count) < 100: continue
		for time in timepoints:
			count = samples.samples_at_time(time)
			if count:
				grav_counts[time][count-1] += 1
				score = (score_motif(pwm, seq1[-20:] + seq2[:8]) - p_min) / (p_max - p_min)

				grav_scores += [score]
		grav_times[samples.num_times()] += 1
		if samples.num_times() == 0: print key

	label = "{}, N = {}".format('Total', sum(total_times[1:]))
	plot_set(total_counts, total_times, label, timepoints, 'b')

	label = "{}, N = {}".format('Weak Motif', sum(i - j for i, j in zip(weak_motif_times[1:], grav_times[1:])))
	plot_set(weak_motif_counts, weak_motif_times, label, timepoints, 'm')

	label = "{}, N = {}".format('Strong Motif', sum(i - j for i, j in zip(motif_times[1:], grav_times[1:])))
	plot_set(motif_counts, motif_times, label, timepoints, 'r')

	label = "{}, N = {}".format('Graveley', sum(grav_times[1:]))
	plot_set(grav_counts, grav_times, label, timepoints, 'g')


	plt.title("Putative Recursive Sites from SJRs")
	plt.ylabel("Number of Recursive Sites")
	plt.xlabel( "    5 min                 10 min                 20 min              total          Number of Timepoints")
	plt.legend(loc=2)

	if directory:
		plt.savefig(directory+'SJRs.png')
	plt.show()

	plt.title('Motif Strength')
	plt.xlabel("    5 min                 10 min                 20 min              total          Number of Timepoints")

	x_vals = []
	labels  = []
	scores = []
	for i, time in enumerate(timepoints):
		begin = i * 4 + 1
		labels += map(str, range(1, len(novel_scores[time]) + 1))
		x_vals += range(begin, begin + len(novel_scores[time]))
		scores += novel_scores[time]

	begin = (len(timepoints) + 1) * 3 + 1
	x_vals += range(begin, begin + len(novel_scores_times))
	labels += map(str, range(1, 5))
	scores += novel_scores_times

	begin = (len(timepoints) + 1) * 3 + 2 + len(novel_scores_times)
	x_vals += [begin, begin + 2]
	labels += ['grav', 'random']
	scores += [grav_scores, random_scores]

	plt.boxplot(scores, positions = x_vals, labels = labels)
	if directory:
		plt.savefig(directory+'numsamples_vs_motif.png')
	plt.show()

	return grav_times, total_times

if __name__ == "__main__":
	import sys
	from intron_rs import IntronRS
	data = open(sys.argv[1], 'r')
	rs = [IntronRS(line) for line in data]
	print main(rs, True)
