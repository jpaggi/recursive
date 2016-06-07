from math import sqrt

class Intron:
	def __init__(self, chrom, start, end, strand, expression, coverage = None):
		self.chrom = chrom
		self.strand = strand
		self.start = start
		self.end = end
		self.rs = []
		self.exon = []
		self.expression = expression
		self.coverage = coverage

	def compatible(self, chrom, strand, rs):
		return (self.chrom == chrom
			and self.strand == strand
			and self.start < rs < self.end)

	def add(self, rs, exon = False):
		self.rs += [rs]
		self.exon += [exon]

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

	def __str__(self):
		return '\t'.join(map(str, [
			self.chrom,
			self.start,
			self.end,
			'.',
			','.join(map(str, self.rs)),
			self.strand]))
