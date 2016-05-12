
class PeakCall:
	def __init__(self, line):
		a = line.strip().split('\t')
		self.chrom = a[0]
		self.rs = int(a[1])
		self.count = int(a[3])
		self.score = float(a[4])
		self.strand = a[5]

	def same(self, other):
		return (self.chrom == other.chrom
			and self.strand == other.strand
			and self.rs == other.rs)

	def merge(self, other):
		self.score = max(self.score, other.score)

	def __str__(self):
		return '\t'.join(map(str, [self.chrom, self.rs, self.rs+1, self.count, self.score, self.strand]))

import sys
cur = PeakCall(sys.stdin.readline())
for line in sys.stdin:
	call = PeakCall(line)
	if cur.same(call):
		cur.merge(call)
	else:
		print cur
		cur = call
print cur