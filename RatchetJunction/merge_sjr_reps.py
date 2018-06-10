import sys

class RS:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom  = a[0]
		self.start  = int(a[1])
		self.end    = int(a[2])
		sample      = a[3]
		count       = int(a[4])
		self.strand = a[5]

		self.counts = {}
		self.counts[sample] = count

	def same(self, other):
		return (self.chrom  == other.chrom
			and self.strand == other.strand
			and self.start  == other.start
			and self.end    == other.end)

	def merge(self, other):
		self.counts = {i+j for i, j in zip(self.counts, other.counts)}

	def __str__(self):
		return '\t'.join(map(str, [
			self.chrom, 
			self.start,
			self.end,
			'.',
			','.join("{}:{}".format(k, v) for k, v in self.counts.items()),
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
