import matplotlib.pyplot as plt
from core_pipeline.get_motifs import *
from core_pipeline.load_genome import *

grav_file = open('../data/graveley.bed', 'r')
data = open('../data/all_sjr_long.bed', 'r')

random_tpss = map(float, open('../data/random_tpss.csv', 'r').read().split(',')[:-1])

genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))

fp_pwm, tp_pwm, fscores, tscores = make_pwm('../data/anno.ss', genome_seq, scores = True)

fp_min, fp_max = get_min_score(fp_pwm), get_max_score(fp_pwm)
tp_min, tp_max = get_min_score(tp_pwm), get_max_score(tp_pwm)

graveley = {}
for line in grav_file:
	chrom, start, end, name, score, strand = line.strip().split('\t')
	rs = int(start) if strand == '-' else int(end)
	graveley[(chrom, strand, rs)] = (set(), 0, '', '')


novel = {}
for line in data:
	chrom, start, end, samples, sjr, strand, seq1, seq2 = line.strip().split('\t')

	if seq1[-2:] != 'AG': continue


	rs = int(start) if strand == '-' else int(end) + 1

	if (chrom, strand, rs) in graveley:
		sample, total, a, b = graveley[(chrom, strand, rs)]
		map(sample.add, samples.split(','))
		graveley[(chrom, strand, rs)] = (sample, total + int(sjr), seq1, seq2)

	if (chrom, strand, rs) in novel:
		sample, total, seq1, seq2 = novel[(chrom, strand, rs)]
		map(sample.add, samples.split(','))
		novel[(chrom, strand, rs)] = (sample, total + int(sjr), seq1, seq2)
	else:
		sample = set()
		map(sample.add, samples.split(','))
		novel[(chrom, strand, rs)] = (sample, int(sjr), seq1, seq2)


totals = [0] * 12
scores = [[] for i in range(12)]
n = 0
for key in novel:
	samples, total, seq1, seq2 = novel[key]
	totals[len(samples)] += 1
	tps = (score_motif(tp_pwm, seq1[-20:]) - tp_min) / (tp_max - tp_min)
	if tps > .8: 
		n +=1
		print seq1[5:] + seq2[:10]
	scores[len(samples)] += [tps]

scores += [[]]
grav = [0] * 12
for key in graveley:
	samples, total, seq1, seq2 = graveley[key]
	grav[len(samples)] += 1
	if total: scores[-1] += [(score_motif(tp_pwm, seq1[-20:]) - tp_min) / (tp_max - tp_min)]


print grav
print totals
print n

plt.boxplot(scores + [tscores, random_tpss], labels = map(str, range(12)) + ['grav', 'anno', 'random'])
plt.show()