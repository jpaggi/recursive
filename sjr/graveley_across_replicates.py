import matplotlib.pyplot as plt
from core_pipeline.get_motifs import *
from core_pipeline.load_genome import *

grav_file = open('../data/graveley.bed', 'r')
data = open('../data/test_adel/merged.bed', 'r')

random_scores = map(float, open('../data/full_motif_AGGT.csv', 'r').read().split(',')[:-1])

genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))

fp_pwm, tp_pwm= make_pwm('../data/anno.ss', genome_seq)

pwm = tp_pwm + fp_pwm

p_min, p_max = get_min_score(pwm), get_max_score(pwm)

graveley = {}
for line in grav_file:
	chrom, start, end, name, score, strand = line.strip().split('\t')
	rs = int(start) if strand == '-' else int(end)
	graveley[(chrom, strand, rs)] = (set(), 0, '', '')


novel = {}
for line in data:
	chrom, start, end, samples, sjr, strand, seq1, seq2 = line.strip().split('\t')[:8]

	if seq1[-2:] != 'AG' or seq2[:2] != 'GT': continue


	rs = int(start) if strand == '-' else int(end) + 1

	if (chrom, strand, rs) in graveley:
		sample, total, a, b = graveley[(chrom, strand, rs)]
		map(sample.add, samples.split(','))
		graveley[(chrom, strand, rs)] = (sample, total + int(sjr), seq1, seq2)
		print line

	if (chrom, strand, rs) in novel:
		sample, total, seq1, seq2 = novel[(chrom, strand, rs)]
		map(sample.add, samples.split(','))
		novel[(chrom, strand, rs)] = (sample, total + int(sjr), seq1, seq2)
	else:
		sample = set()
		map(sample.add, samples.split(','))
		novel[(chrom, strand, rs)] = (sample, int(sjr), seq1, seq2)


totals = [0] * 13
gt = [0] * 13
scores = [[] for i in range(13)]
n = 0
for key in novel:
	samples, total, seq1, seq2 = novel[key]
	totals[len(samples)] += 1
	tps = (score_motif(pwm, seq1[-20:] + seq2[:8]) - p_min) / (p_max - p_min)
	if tps > .875:
		if key not in graveley:
			n +=1
	if seq2[:2] == 'GT':
		scores[len(samples)] += [tps]
		gt[len(samples)] += 1


scores += [[]]
grav = [0] * 13
for key in graveley:
	samples, total, seq1, seq2 = graveley[key]
	grav[len(samples)] += 1
	if total:
		#print seq1[5:] + seq2[:10]
		scores[-1] += [(score_motif(pwm, seq1[-20:] + seq2[:8]) - p_min) / (p_max - p_min)]


print grav
print totals
print n

plt.title('Number of Samples Per Putative Recursive Site')
plt.xlabel('Number of Samples')
plt.ylabel('Number of Putative Recursive Sites')
plt.xlim(-1, 13)
plt.bar(range(13), totals, align = 'center')
plt.bar(range(13), gt, align = 'center')
plt.show()

plt.title('Motif Strength')
plt.xlabel('Number of Samples')
plt.boxplot(scores + [random_scores], labels = map(str, range(13)) + ['grav', 'random'])
plt.show()