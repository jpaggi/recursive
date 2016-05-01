from load_genome import revcomp
from math import log

FP_LEN = 8
TP_LEN = 20

def get_pwm(seqs):
	length = len(seqs[0])
	pwm = [{'A': 1, 'C': 1, 'G': 1, 'T': 1} for i in range(length)]
	for seq in seqs:
		for i, char in enumerate(seq):
			if char not in pwm[i]: continue
			pwm[i][char] += 1

	for pos in range(length):
		for char in pwm[pos]:
			pwm[pos][char] = pwm[pos][char] / float(len(seqs))

	return pwm

def get_scores(seqs, pwm):
	min_score, max_score = get_min_score(pwm), get_max_score(pwm)
	scores = []
	for seq in seqs:
		scores += [(score_motif(pwm, seq) - min_score) / (max_score - min_score)]
	return scores


def make_pwm(ss, genome, scores = False):
	ss = open(ss, 'r')

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
			fps += [genome[chrom][pos+1:pos+FP_LEN+1]]

	for chrom in threes['+']:
		for pos in threes['+'][chrom]:
			tps += [genome[chrom][pos-TP_LEN: pos]]


	for chrom in fives['-']:
		for pos in fives['-'][chrom]:
			fps += [revcomp(genome[chrom][pos-FP_LEN: pos])]

	for chrom in threes['-']:
		for pos in threes['-'][chrom]:
			tps += [revcomp(genome[chrom][pos+1:pos+TP_LEN+1])]


	fpwm, tpwm = get_pwm(fps), get_pwm(tps)
	if scores:
		fscores, tscores = get_scores(fps, fpwm), get_scores(tps, tpwm)
		return fpwm, tpwm, fscores, tscores
	return fpwm, tpwm

def get_min_score(pwm):
	bits = 0
	for w in pwm:
		p = min([w[c] for c in w])
		bits += log(p / 0.25, 2)
	return bits

def get_max_score(pwm):
	bits = 0
	for w in pwm:
		p = max([w[c] for c in w])
		bits += log(p / 0.25, 2)
	return bits


def score_motif(pwm, seq):
	assert len(pwm) == len(seq)
	bits = 0
	for w, c in zip(pwm, seq):
		if c == 'N': return 0
		bits += log(w[c] / 0.25, 2)
	return bits

def search_for_motif(pwm, seq, min_score, max_score, AG = False):
	best = (-1, 0)
	for i in range(len(seq) - len(pwm) + 1):
		if AG and seq[i + TP_LEN - 2: i + TP_LEN] != 'AG': continue
		p = score_motif(pwm, seq[i:i+len(pwm)])
		if p > best[1]:
			best = (i, p)

	return best[0], (best[1] - min_score) / (max_score - min_score)


