import sys

SAMPLES = ['sjr_5', 'pe_5', 'sjr_10', 'pe_10', 'sjr_20', 'pe_20', 'sjr_total', 'pe_total']
index = lambda x: SAMPLES.index(x)

class RS:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom  = a[0]
		self.start  = int(a[1])
		self.end    = int(a[2])
		sample      = a[3][:-5]
		count       = int(a[4])
		self.strand = a[5]

		self.counts = [0] * len(SAMPLES)
		self.counts[index(sample)] += count

	def same(self, other):
		return (self.chrom  == other.chrom
			and self.strand == other.strand
			and self.start  == other.start
			and self.end    == other.end)

	def merge(self, other):
		self.counts = [i+j for i, j in zip(self.counts, other.counts)]

	def __str__(self):
		return '\t'.join(map(str, [
			self.chrom, 
			self.start,
			self.end,
			'.',
			','.join(map(str, self.counts)),
			self.strand]))

cur = RS(sys.stdin.readline())
for line in sys.stdin:
	rs = RS(line)
	if cur.same(rs):
		cur.merge(rs)
	else:
		print cur
		cur = rs
print cur
