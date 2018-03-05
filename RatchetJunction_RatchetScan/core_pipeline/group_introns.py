import sys

if sys.argv[1] == '-':
	reads = sys.stdin
else:
	reads = open(sys.argv[1], 'r')

SAMPLE = sys.argv[2]

class Group:
	def __init__(self, read):
		a = read.strip().split('\t')
		self.chrom = a[0]
		self.start = int(a[1])
		self.end = int(a[2])
		self.name = a[3]
		blocks = eval(a[4])
		self.offsets = set()
		self.offsets.add((blocks[0][0], blocks[-1][1]))
		self.strand = a[5]
		self.total = 1

	def same(self, other):
		return (self.chrom  == other.chrom
			and self.start  == other.start
			and self.end    == other.end
			and self.strand == other.strand)

	def merge(self, other):
		self.offsets |= other.offsets
		self.total += other.total

	def __str__(self):
		return '\t'.join(map(str, [self.chrom, self.start, self.end, SAMPLE, self.total, self.strand]))

cur = Group(reads.readline())
for read in reads:
	entry = Group(read)
	if cur.same(entry):
		cur.merge(entry)
	else:
		print cur
		cur = entry
print cur
