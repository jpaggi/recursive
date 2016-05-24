from math import sqrt

class Intron:
	def __init__(self, chrom, start, end, strand, expression):
		self.chrom = chrom
		self.strand = strand
		self.start = start
		self.end = end
		self.rs = []
		self.expression = expression

	def compatible(self, chrom, strand, rs):
		return (self.chrom == chrom
			and self.strand == strand
			and self.start < rs < self.end)

	def add(self, rs):
		self.rs += [rs]

	def num_rs(self):
		return len(self.rs)

	def length(self):
		return self.end - self.start

	def rs_lengths(self):
		sites = sorted(self.rs)
		sites = [self.start] + sites + [self.end]
		lengths = [e - b for b, e, in zip(sites[:-1], sites[1:])]
		return [e - b for b, e, in zip(sites[:-1], sites[1:])]

	def rs_length_dispersion(self):
		lengths = self.rs_lengths()
		expectation = lambda x: sum(x) / float(len(x))
		average = expectation(lengths)
		return sqrt(expectation(map(lambda x: (x - average) ** 2 , lengths))) / average
