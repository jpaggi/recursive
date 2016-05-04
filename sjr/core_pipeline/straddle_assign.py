from load_genome import load_genome, revcomp
from get_motifs import make_pwm, score_motif, get_min_score, get_max_score, search_for_motif
import sys

genome_seq = load_genome(open(sys.argv[4], 'r'))
fp_pwm, tp_pwm = make_pwm(sys.argv[3], genome_seq)

pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

sjr_file = open(sys.argv[2])

class Group:
	def __init__(self, chrom, five, three, strand, sample, five_insert = None, count = None):
		assert (five_insert == None) != (count == None)
		self.strand = strand
		self.chrom = chrom
		self.five = five
		self.sample = sample

		if five_insert != None:
			self.range = self._get_range(three, five_insert)
			self.sjr = 0
			self.per = 1
		else:
			self.sjr = count
			self.per = 0
			self.range = three

	def _get_range(self, three, five_insert):
		if strand == '+':
			start = three + five_insert - 300
			end = three
		else:
			start = three
			end = three + 300 - five_insert
		return (start, end)

	def compatible(self, three, five_insert):
		r = self._get_range(three, five_insert)
		if self.sjr and r[0] <= self.range <= r[1]:
			return self.sjr + self.per
		elif not self.sjr and r[0] <= self.range[1] and r[1] >= self.range[0]:
			return self.sjr + self.per
		else:
			return -1

	def add(self, three, five_insert):
		assert self.compatible(three, five_insert)
		self.per += 1
		if not self.sjr:
			r = self._get_range(three, five_insert)
			self.range = (max(r[0], self.range[0]), min(r[1], self.range[1]))

	def _get_ss(self):
		if self.sjr:
			return self.range
		elif self.strand ==  '+':
			start = self.range[0] - len(tp_pwm) - 50
			end = self.range[1] + len(fp_pwm) + 10
			seq = genome_seq[self.chrom][start:end]
			i, p = search_for_motif(pwm, seq, min_score, max_score, AG = True)
			if i == -1: return None
			return start + i + len(tp_pwm)
		else:
			start = self.range[0] - len(fp_pwm) - 10
			end = self.range[1] + len(tp_pwm) + 50
			seq = revcomp(genome_seq[self.chrom][start:end])
			i, p = search_for_motif(pwm, seq, min_score, max_score, AG = True)
			if i == -1: return None
			return end - i - len(tp_pwm)

	def print_entry(self):
		rs = self._get_ss()
		if self.strand == '+':
			start, end = self.five, rs 
		else:
			start, end = rs, self.five
		if rs and self.per:
			print '\t'.join(map(str, [self.chrom, start, end, self.sample, self.per, self.strand, self.range]))


sjr = {}
for line in sjr_file:
	chrom, start, end, samples, count, strand = line.strip().split('\t')[:6]

	if strand == '+':
		five, three = int(start), int(end)
	else:
		five, three = int(end), int(start)

	key = (chrom, five, strand)

	if key in sjr:
		sjr[key] += [(three, int(count))]
	else:
		sjr[key] = [(three, int(count))]

per_file = open(sys.argv[1], 'r')

per = {}

for line in per_file:
	chrom, inner_left, inner_right, sample, five, strand = line.strip().split('\t')
	inner_left, inner_right, five = int(inner_left), int(inner_right), int(five)


	if strand ==  '+':
		five_insert = five - inner_left
		three = inner_right
	else:
		five_insert = inner_right - five
		three = inner_left

	key = (chrom, five, strand)

	if key in per:
		per[key] += [(three, five_insert)]
	else:
		per[key] = [(three, five_insert)]


for key in per:
	chrom, five, strand = key
	groups = [Group(chrom, five, three, strand, sample, count = count) for three, count in sjr[key]] if key in sjr else []
	for three, five_insert in per[key]:
		best = (None, 0)
		for i, group in enumerate(groups):
			count = group.compatible(three, five_insert)
			if count > best[0]:
				best = (i, count)
		if best[0]:
			groups[best[0]].add(three, five_insert)
		else:
			groups += [Group(chrom, five, three, strand, sample, five_insert = five_insert)]
	for group in groups: 
		group.print_entry()
