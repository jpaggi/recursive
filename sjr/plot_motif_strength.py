import sys

graveley_file = open('../data/recursiveintrons.bed')
from load_genome import load_genome, revcomp
from get_motifs import make_pwm, score_motif, get_min_score, get_max_score
import sys
import matplotlib.pyplot as plt
THRESH = 0.87

genome_seq = load_genome(open(sys.argv[1], 'r'))
fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)

fp_min, fp_max = get_min_score(fp_pwm), get_max_score(fp_pwm)
tp_min, tp_max = get_min_score(tp_pwm), get_max_score(tp_pwm)

graveley = []
for line in graveley_file:
	chrom, start, end, igv, rs, strand = line.strip().split('\t')
	rs = map(int, rs.split(','))
	start, end = int(start), int(end)
	graveley += [(chrom, start, end, strand, rs)]


total = [0] * 9
n_above_thresh = [0] * 9
g_above_thresh = [0] * 9
min_g = 1
for line in sys.stdin:
	chrom, start, end, reps, offsets, strand, seq1, seq2 = line.strip().split('\t')
	if seq1[-2:] != 'AG': continue
	start, end = int(start), int(end)
	reps = reps.split(',')
	reps = map(lambda x: x.split('_rep'), reps)
	offsets = int(offsets)
	pwm = tp_pwm + fp_pwm

	min_score = fp_min + tp_min
	max_score = fp_max + tp_max
	seq = seq1[-15:] + seq2[:8]
	score = (score_motif(pwm, seq) - min_score) / (max_score - min_score)

	five = 0
	ten = 0
	twenty = 0
	for rep in reps:
		if rep[0] == '5':
			five += 1
		elif rep[0] == '10':
			ten += 1
		else:
			twenty += 1

	seen = False
	for (g_chrom, g_start, g_end, g_strand, g_rs) in graveley:
		# - strand indexing the same
		# + strand g_rs == end + 1
		if g_chrom == chrom and g_strand == strand:# and #(g_start <= start <= g_end or g_start <= end <= g_end):
			if strand == '+':
				seen = seen or (end + 1 in g_rs)
				rs = end
			else:
				seen = seen or (start in g_rs)
				rs = start

	total[five + ten + twenty - 1] += 1

	if not seen:
		#print chrom, start, end, strand, seq1, seq2, rs, 'novel'
		plt.scatter(five, score, c = 'r')
		if score > THRESH:
			n_above_thresh[five + ten + twenty - 1] += 1
	else:
		#print chrom, start, end, strand, seq1, seq2, rs, 'graveley'
		plt.scatter(five, score, c = 'g')
		min_g = min(min_g, score)
		if score > THRESH:
			g_above_thresh[five + ten + twenty - 1] += 1

print 'threshold = ' + str(THRESH)
print 'novel', n_above_thresh, sum(n_above_thresh)
print 'gravely', g_above_thresh, sum(g_above_thresh)
print 'total', total, sum(total) 
print min_g
plt.show()
