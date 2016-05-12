import sys

class IntronRS:
	def __init__(self, line):
		attributes = line.strip().split('\t')

		self.chrom   = attributes[0]
		self.start   = int(attributes[1])
		self.end     = int(attributes[2])
		self.samples = set(attributes[3].split(','))
		self.rs      = int(attributes[4])
		self.strand  = attributes[5]     
		self.sjr     = int(attributes[6])
		self.intron  = int(attributes[7])
		self.seq1    = attributes[8]
		self.seq2    = attributes[9]
		self.score1  = float(attributes[10])
		self.score2  = float(attributes[11])

	def same(self, other):
		return (self.chrom  == other.chrom
			and self.start  == other.start
			and self.end    == other.end
			and self.strand == other.strand
			and self.rs     == other.rs)

	def merge(self, other):
		assert self.same(other)
		self.sjr += other.sjr
		self.samples |= other.samples

	def __str__(self):
		return '\t'.join(map(str, [
				self.chrom, 
				self.start,
				self.end,
				','.join(self.samples),
				self.rs,
				self.strand,
				self.sjr,
				self.intron,
				self.seq1,
				self.seq2,
				self.score1,
				self.score2
			]))



cur = IntronRS(sys.stdin.readline())
for line in sys.stdin:
	intron = IntronRS(line)
	if cur.same(intron):
		cur.merge(intron)
	else:
		print cur
		cur = intron
print cur