import sys
from load_genome import *
from get_motifs import *


SAMPLES = ["{}_rep{}".format(a, b) for a in ['5', '10', '20'] for b in ['A', 'B', 'C']]
SAMPLES += ['total_repA', 'total_repB']

def MAP(x):
	if x == '5': return 0
	elif x == '1': return 1
	elif x == '2': return 2
	elif x == 't': return 3
	assert False


files = {}
for i in SAMPLES:
	files[i] = open("non_seq_data/rec_to_rec_{}.bed".format(i))


sites = {}

for file in files:
	for line in files[file]:
		chrom, start, end, ID, coords, strand = line.strip().split('\t')
		start, end = int(start), int(end)

		key = (chrom, start, end, strand)

		if not key in sites:
			sites[key] = [set(), set(), set(), set()]
		sites[key][MAP(file[0])].add(coords)


genome_seq = load_genome(open('annotations/dmel/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
fp_pwm, tp_pwm = make_pwm('annotations/dmel/dmel-anno-splice-sites.ss', genome_seq)

pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)


for site in sites:
	chrom, start, end, strand = site

	sites[site] = map(len, sites[site])
	sjr = ','.join(map(str, sites[site]))


	seq1 = genome_seq[chrom][start-30:start+30]
	seq2 = genome_seq[chrom][end-30:end+30]
	if strand == '-':
		seq1 = revcomp(seq1)[10:38]
		seq2 = revcomp(seq2)[11:39]
	else:
		seq1 = seq1[9:37]
		seq2 = seq2[10:38]

	score1 = (score_motif(pwm, seq1) - min_score) / (max_score - min_score)
	score2 = (score_motif(pwm, seq2) - min_score) / (max_score - min_score)

	#if strand == '+' and score1 < .85: continue
	#if strand == '-' and score2 < .85: continue

	if score1 < .85 or score2 < .85: continue

	if sum(sites[site]) < 3: continue

	print '\t'.join(map(str, [chrom, start, end, score1, score2, strand, sjr]))


