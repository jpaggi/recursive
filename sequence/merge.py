import sys

class ExpBed:
	def __init__(self, line):
		chrom, start, end, name, sjr, strand, expression = line.strip().split('\t')
		self.chrom = chrom
		self.start = int(start)
		self.end = int(end)
		self.name = name
		self.sjr = int(sjr)
		self.strand = strand
		self.expression = map(int, expression.split(','))

	def same(self, other):
		return (
			self.chrom  == other.chrom and
			self.start  == other.start and
			self.end    == other.end   and
			self.strand == other.strand
			)

	def add(self, other):
		assert self.chrom == other.chrom
		assert self.start == other.start
		assert self.end == other.end
		assert self.strand == other.strand
		assert len(self.expression) == len(other.expression)
		self.sjr += other.sjr
		self.expression = [s+o for s, o in zip(self.expression, other.expression)]

	def __str__(self):
		express = ','.join(map(str, self.expression))
		return '\t'.join(map(str, [self.chrom, self.start, self.end, self.name, self.sjr, self.strand, express]))

cur_entry = ExpBed(sys.stdin.readline())
for line in sys.stdin:
	entry = ExpBed(line)

	if cur_entry.same(entry):
		cur_entry.add(entry)
	else:
		print cur_entry
		cur_entry = entry
print cur_entry
