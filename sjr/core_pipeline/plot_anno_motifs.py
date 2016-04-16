import sys
#import matplotlib.pyplot as plt
from get_motifs import make_pwm, get_min_score, get_max_score, score_motif
from load_genome import load_genome, revcomp

print 'currently making scrambled motifs'

ss = open(sys.argv[1], 'r')
genome = load_genome(open(sys.argv[2], 'r'))

fives = {'+': {}, '-': {}}
threes = {'+': {}, '-': {}}

for line in ss:
	chrom, start, end, strand = line.strip().split('\t')
	chrom = chrom[3:]
	if strand == '+':
		five, three = int(start), int(end)
	else:
		three, five = int(start), int(end)

	if chrom in fives[strand]:
		fives[strand][chrom].add(five)
	else:
		fives[strand][chrom] = set([five])

	if chrom in threes[strand]:
		threes[strand][chrom].add(three)
	else:
		threes[strand][chrom] = set([five])

fps = []
tps = []
for chrom in fives['+']:
	for pos in fives['+'][chrom]:
		pos += 100	
		fps += [genome[chrom][pos+1:pos+9]]

for chrom in threes['+']:
	for pos in threes['+'][chrom]:
		pos -= 100
		tps += [genome[chrom][pos-15: pos]]


for chrom in fives['-']:
	for pos in fives['-'][chrom]:
		pos -= 100
		fps += [revcomp(genome[chrom][pos-8: pos])]

for chrom in threes['-']:
	for pos in threes['-'][chrom]:
		pos += 100
		tps += [revcomp(genome[chrom][pos+1:pos+16])]




fp_pwm, tp_pwm = make_pwm(sys.argv[1], genome)

fp_min, fp_max = get_min_score(fp_pwm), get_max_score(fp_pwm)
tp_min, tp_max = get_min_score(tp_pwm), get_max_score(tp_pwm)

fps_scores = []

for seq in fps:
	fps_scores += [(score_motif(fp_pwm, seq[:8]) - fp_min) / (fp_max - fp_min)]

tps_scores = []

for seq in tps:
	tps_scores += [(score_motif(tp_pwm, seq[-15:]) - tp_min) / (tp_max - tp_min)]

f_c = 0
for score in fps_scores:
	if score > .8:
		f_c += 1
t_c = 0
for score in tps_scores:
	if score > .8:
		t_c += 1

print f_c / float(len(fps_scores)), t_c / float(len(tps_scores))


#plt.hist(fps_scores, normed=True, cumulative=True, bins = 40)
#plt.show()
