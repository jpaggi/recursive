import sys
from get_motifs import *
from load_genome import *

class Site:
	def __init__(self, chrom, strand, rs):
		self.chrom = chrom
		self.strand = strand
		self.rs = rs
		self.motifscore = self.set_motif_score()

		self.sjr = {}
		self.saw_prob = 0
		self.saw_score = 0

	def assign_intron(self, five, three, expression):
		self.five = five
		self.three = three
		self.intron_expression = expression

	def set_motif_score(self):
		if self.strand == '+':
			motif = genome_seq[self.chrom][self.rs - len(tp_pwm): self.rs + len(fp_pwm)]
		else:
			motif = genome_seq[self.chrom][self.rs - len(fp_pwm): self.rs + len(tp_pwm)]
			motif = revcomp(motif)
		return (score_motif(pwm, motif) - min_score) / (max_score - min_score)

	def add_sjr(self, five, counts):
		for sample in counts.split(','):
			name, count = sample.split(':')
			if five not in self.sjr: self.sjr[five] = {}
			assert name not in self.sjr[five]
			self.sjr[five][name] = int(count)

	def sjr_str(self):
		return ' '.join(str(five) + ':' + ','.join("{}:{}".format(name, count) for name, count in counts.items())
						 for five, counts in self.sjr.items()
						 )

	def add_sawtooth(self, prob, score):
		self.saw_prob = prob
		self.saw_score = score

	def __str__(self):
		return '\t'.join(map(str, [self.chrom, self.rs, self.motifscore,
			                       self.sjr_str(), self.saw_prob, self.saw_score]))


genome_seq = load_genome(open(sys.argv[4]))
fp_pwm, tp_pwm = make_pwm(sys.argv[5], genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)
sites = {}

# RatchetJunction
with open(sys.argv[1]) as fp:
	for line in fp:
		chrom, start, end, _, counts, strand = line.strip().split()
		if strand == '+':
			five, rs = int(start), int(end)
		else:
			five, rs = int(end), int(start)
		key = (chrom, strand, rs)

		if key not in sites: sites[key] = Site(chrom, strand, rs)
		sites[key].add_sjr(five, counts)

# RatchetPair
with open(sys.argv[2]) as fp:
	for line in fp:
		chrom, start, end, pi, counts, strand = line.strip().split()
		start, end = int(start), int(end)
		if strand == '+':
			five, rs = int(start), int(end)
		else:
			five, rs = int(end), int(start)
		key = (chrom, strand, rs)

		if key not in sites: sites[key] = Site(chrom, strand, rs)
		sites[key].add_sjr(five, counts)

# RatchetScan
with open(sys.argv[3]) as fp:
	for line in fp:
		chrom, start, end, prob, score, strand = line.strip().split()
		prob, score = float(prob), float(score)
		key = (chrom, strand, int(start))

		if key not in sites: sites[key] = Site(chrom, strand, rs)
		sites[key].add_sawtooth(prob, score)

for site in sites.values():
	print site
