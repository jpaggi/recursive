from math import sqrt

class Intron:
	def __init__(self, intronrs):
		self.rs = [intronrs]
		self.expression = intronrs.intron
		self.start = intronrs.start
		self.end = intronrs.end
		self.strand = intronrs.strand
		self.chrom = intronrs.chrom

	def five(self):
		return self.start if self.strand == '+' else self.end

	def three(self):
		return self.end if self.strand == '+' else self.start

	def compatible(self, other):
		return (self.chrom == other.chrom
			and self.strand == other.strand
			and self.five() == other.five())

	def add(self, other):
		assert self.compatible(other)
		for rs in other.rs:
			if not rs.rs in map(lambda x: x.rs, self.rs):
				self.rs += other.rs

		if self.strand == '+':
			self.end = max(self.three(), other.three())
		else:
			self.start = min(self.three(), other.three())

	def num_rs(self):
		return len(self.rs)

	def length(self):
		return self.end - self.start

	def rs_lengths(self):
		sites = sorted(map(lambda x: x.rs, self.rs))
		sites = [self.start] + sites + [self.end]
		lengths = [e - b for b, e, in zip(sites[:-1], sites[1:])]
		return [e - b for b, e, in zip(sites[:-1], sites[1:])]

	def rs_length_dispersion(self):
		lengths = self.rs_lengths()
		expectation = lambda x: sum(x) / float(len(x))
		average = expectation(lengths)
		return sqrt(expectation(map(lambda x: (x - average) ** 2 , lengths))) / average
