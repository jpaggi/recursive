from load_genome import load_genome, revcomp
from get_motifs import make_pwm, score_motif, get_min_score, get_max_score, search_for_motif
import sys
from math import exp
from intron_rs import IntronRS

"""
1) SJR file
2) PER file
3) motif junction file
4) genome fasta file
"""

SAMPLES = ['sjr_5', 'pe_5', 'sjr_10', 'pe_10', 'sjr_20', 'pe_20', 'sjr_total', 'pe_total']
index = lambda x: SAMPLES.index(x)

genome_seq = load_genome(open(sys.argv[3], 'r'))
fp_pwm, tp_pwm = make_pwm(sys.argv[2], genome_seq)

pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

class Intron:
	def __init__(self, chrom, five, strand):
		self.chrom = chrom
		self.five = int(five)
		self.strand = strand

        # things that will be altered by EM
		self.per = []
		self.sjr = {}
		self.pi = None
		self.seq = None
		self.a_m = None
		self.scores = None
		self.counts = None
		self.assignments = None
		self.f_a_s = lambda x: 0 if sum(x) < 2 else sum(x) / 20.0
		self.a_s = None

		# fixed parameters
		self.mu = 1.0
		self.insert_u = 120.0
		self.insert_o = 80.0
		self.motif_cutoff = 0.83
		self.convergence = 0.0000000001
		self.max_insert = int(self.insert_o * 15 + self.insert_u)
		self.allowed_overhang = 15

	def key(self):
		return (self.chrom, self.five, self.strand)

	def add_sjr(self, three, count):
		i = abs(self.five - three)
		self.sjr[i] = count

	def add_per(self, three, five_insert, sample):
		i = abs(self.five - three)
		self.per += [(i, five_insert, sample)]

	def length(self):
		return abs(self.five - self._get_downstream_end())

	def _get_downstream_end(self):
		return max(map(lambda x: x[0], self.per) + self.sjr.keys())

	def _set_seq(self):
		if self.strand == '+':
			start, end = self.five, self.five + self._get_downstream_end() + 10
		else:
			start, end = self.five - self._get_downstream_end() - 10, self.five

		start, end = start - len(tp_pwm), end + len(tp_pwm)
		self.seq = genome_seq[self.chrom][start:end]
		if self.strand == '-': self.seq = revcomp(self.seq)

	def _set_a_m(self):
		best = max(self.scores)
		if best == 0:
			self.a_m = [0] * len(self.scores)
		else:
			self.a_m = []
			for score in self.scores:
				self.a_m += [self.a_s * self.mu * (score - self.motif_cutoff) / best]

	def _set_scores(self):
		self.scores = []
		for i in xrange(len(tp_pwm), len(self.seq) - len(fp_pwm)):
			motif = self.seq[i-len(tp_pwm):i + len(fp_pwm)]
			if motif[len(tp_pwm) - 2: len(tp_pwm) + 2] != 'AGGT':
				self.scores += [0]
			else:
				self.scores += [(score_motif(pwm, motif) - min_score) / (max_score - min_score)]

	def _set_pi(self):
		self.pi = [10 * max(0,  score - self.motif_cutoff) / float(len(self.scores)) for score in self.scores]

	def _e_step(self):
		"""
		Here you update self.counts based on fractional read assignments.
		"""
		self.counts = [0] * len(self.scores)

		# sjr are assigned deterministically
		for pos in self.sjr:
		    self.counts[pos] += self.sjr[pos]

		# per are assigned based on prob of generation
		for three, five_insert, sample in self.per:
		    total = 0.0
		    for pos in xrange(max(0, three - self.max_insert), three + self.allowed_overhang):
		        if self.pi[pos] == 0: continue
		        inferred_insert = five_insert + three - pos
		        total += self._prob_insert(inferred_insert) * self.pi[pos]
		    if total == 0: continue # remove read from consideration.....
		    for pos in xrange(max(0, three - self.max_insert), three + self.allowed_overhang):
		        if self.pi[pos] == 0: continue
		        inferred_insert = five_insert + three - pos
		        self.counts[pos] += (self._prob_insert(inferred_insert) * self.pi[pos]) / total

	def _assign_per(self):
		self.assignments = []
		for three, five_insert, sample in self.per:
			best = (-1, 0)
			for pos in xrange(max(0, three - self.max_insert), three + self.allowed_overhang):
				inferred_insert = five_insert + three - pos
				p = self._prob_insert(inferred_insert) * self.pi[pos]
				if p > best[1]:
					best = (pos, p)
			if best[0] != -1:
				self.assignments += [(three, five_insert, best[0], sample)]

	def _m_step(self):
		pi = [max(0, self.counts[i] - self.a_s + self.a_m[i]) for i in xrange(len(self.counts))]
		total = sum(pi)
		if total == 0:
			return self.pi # signal to end iteration
		return map(lambda x: x / total, pi)

	def _em(self):
		self.a_s = 1 + (len(self.sjr) + len(self.per))
		self._set_seq()
		self._set_scores()
		self._set_pi()

		done = False
		while not done:
			self._e_step()
			self.a_s = self.f_a_s(self.counts)
			self._set_a_m()
			pi = self._m_step()
			done = not any(abs(pi[i]-self.pi[i]) > self.convergence for i in xrange(len(pi)))
			self.pi = pi
		self._assign_per()


	def _prob_insert(self, insert):
		dev = abs(insert - self.insert_u)
		if dev > self.max_insert: return 0
		return exp(dev ** 2 / (-2* self.insert_o ** 2))

	def i_to_g(self, pos):
		if self.strand == '+':
			return self.five + pos
		else:
			return self.five - pos

	def get_assignments(self):
		self._em()
		ends = {}
		for three, insert, prediction, sample in self.assignments:
			if prediction not in ends:
				ends[prediction] = [0] * len(SAMPLES)
			ends[prediction][index(sample)] += 1
		return [(self.chrom, self.five, self.i_to_g(end), ','.join(map(str, ends[end])),
				 self.pi[end], self.strand) for end in ends]


per_file = sys.stdin
introns = {}
for line in per_file:
	chrom, inner_left, inner_right, sample, five, strand = line.strip().split('\t')
	inner_left, inner_right, five = int(inner_left), int(inner_right), int(five)
	if strand ==  '+':
		five_insert = five - inner_left
		three = inner_right
	else:
		five_insert = inner_right - five
		three = inner_left
	sample = '_'.join(sample.split('_')[:2])
	key = (chrom, five, strand)

	if key not in introns: introns[key] = Intron(chrom, five, strand)

	introns[key].add_per(three, five_insert, sample)

sjr_file = open(sys.argv[1])
for line in sjr_file:
	sjr = IntronRS(line, genome_seq)
	if not sjr.aggt(): continue
	if sjr.key() in introns: introns[sjr.key()].add_sjr(sjr.three(), sjr.read_count())

for intron in introns:
	assignments = introns[intron].get_assignments()
	for chrom, five, three, count, pi, strand in assignments:
		if strand == '+': start, end = five, three
		else: start, end = three, five
		print '\t'.join(map(str, [chrom, start, end, pi, count, strand]))
