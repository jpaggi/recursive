

class Read:
	def __init__(self, entry):
		qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = entry.strip().split('\t')[:11]

		self.qname = qname
		self.rname = rname
		self.pos = int(pos) - 1 if pos != '*' else -1
		self.cigar = cigar
		self.seq = seq
		self.qual = qual

		flag = "{0:b}".format(int(flag))

		self.is_reverse = (len(flag) > 4 and flag[-5] == '1')
		self.is_read1 = (len(flag) > 6 and flag[-7] == '1')

	def align_str(self):
		strand = '-' if self.is_reverse else '+'
		return ','.join([self.rname, str(self.pos), self.cigar, strand])


class SamIter:
	def __init__(self, sam, unmapped = False):
		self.sam = open(sam, 'r')
		self.unmapped = unmapped

	def __iter__(self):
		return self

	def next(self):
		line = self.sam.readline()

		if not line:
			raise StopIteration

		while line[0] == '@':
			line = self.sam.readline()

		if line.split('\t')[2] != '*' or self.unmapped:
			return Read(line)
		else:
			return self.next()

