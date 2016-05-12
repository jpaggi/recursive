"""
take as input the loj of
"""
import sys

class Intersect:
	def __init__(self, line):
		a = line.strip().split('\t')

		self.chrom = a[0]
		self.start = int(a[1])
		self.end = int(a[2])
		self.name = a[3]
		self.count = int(a[10])
		self.strand = a[5]

	def same(self, other):
		return (self.chrom == other.chrom
			and self.strand == other.strand
			and self.start == other.start
			and self.end == other.end)

	def merge(self, other):
		self.count += other.count

	def __str__(self):
		return '\t'.join(map(str, 
			[self.chrom,
			self.start,
			self.end,
			self.name,
			self.count,
			self.strand]))

cur = Intersect(sys.stdin.readline())
for line in sys.stdin:
	intersect = Intersect(line)
	if cur.same(intersect):
		cur.merge(intersect)
	else:
		print cur
		cur = intersect
print cur
