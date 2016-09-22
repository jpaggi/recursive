import sys
from load_genome import load_genome, revcomp
from get_motifs import make_pwm, get_min_score, get_max_score, score_motif
#import matplotlib.pyplot as plt

BINS = 100
TP_LEN = 20
AG = True
MIN_INTRON_LEN = 10000


fasta_file = sys.argv[1]
ss_file = sys.argv[2]
bed = open(sys.argv[3], 'r')
out = open(sys.argv[4], 'w')

genome = load_genome(open(fasta_file, 'r'))
fp_pwm, tp_pwm = make_pwm(ss_file, genome)

pwm = tp_pwm + fp_pwm

min_score, max_score = get_min_score(pwm), get_max_score(pwm)
most = 0
scores = [0] * BINS


eighty, eighty_five, ninety = 0, 0, 0

for line in bed:
	chrom, start, end, strand = line.strip().split('\t')
	start, end = int(start), int(end)

	if end - start < MIN_INTRON_LEN: continue

	seq = genome[chrom][start+20: end-20]
	if strand == '-': seq = revcomp(seq)

	total = 0

	for i in range(0, len(seq) - len(pwm)):
		if AG and seq[i+TP_LEN - 2:i+TP_LEN+2] != 'AGGT': continue

		score = (score_motif(pwm, seq[i:i+len(pwm)]) - min_score) / float(max_score - min_score)

		scores[int(score * BINS)] += 1

		pos = start + 20 + TP_LEN + i if strand == '+' else end - 20 - TP_LEN - i
		out.write('\t'.join(map(str, [chrom, start, end, pos, score, strand]))+ '\n' )

		total += score > .8

		eighty += score > .8
		eighty_five += score > .85
		ninety += score > .9

		#print genome[chrom][pos-2:pos+2]

	most = max(most, total)

print eighty, eighty_five, ninety
print most
