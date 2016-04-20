import sys
from load_genome import load_genome, revcomp
from get_motifs import make_pwm, get_min_score, get_max_score, score_motif
import matplotlib.pyplot as plt

BINS = 100
TP_LEN = 20
AG = False
MIN_INTRON_LEN = 10000


fasta_file = sys.argv[1]
ss_file = sys.argv[2]
bed = open(sys.argv[3], 'r')

genome = load_genome(open(fasta_file, 'r'))
fp_pwm, tp_pwm = make_pwm(ss_file, genome)

pwm = tp_pwm + fp_pwm

min_score, max_score = get_min_score(pwm), get_max_score(pwm)

scores = [0] * BINS
for line in bed:
	chrom, start, end, name, name2, strand = line.strip().split('\t')
	chrom, start, end = chrom[3:], int(start), int(end)

	if end - start < MIN_INTRON_LEN: continue

	seq = genome[chrom][start+20: end-20]
	if strand == '-': seq = revcomp(seq)

	for i in range(0, len(seq) - len(pwm), 5):
		if AG and seq[i+TP_LEN - 2:i+TP_LEN] != 'AG': continue

		score = (score_motif(pwm, seq[i:i+len(pwm)]) - min_score) / float(max_score - min_score)

		scores[int(score * BINS)] += 1

total = float(sum(scores))
print total
print scores
scores = map(lambda x: x / total, scores)

print scores

plt.plot(scores)
plt.show()

