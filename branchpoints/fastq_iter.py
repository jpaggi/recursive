def complement(char):
	if char == 'A': return 'T'
	if char == 'C': return 'G'
	if char == 'G': return 'C'
	if char == 'T': return 'A'
	return char

def revcomp(seq):
	return ''.join(map(complement, seq[::-1]))


class FastQIter:
	def __init__(self, reads1, reads2):
		self.reads1 = open(reads1, 'r')
		self.reads2 = open(reads2, 'r')
		self.last_ID1 = ''
		self.last_ID2 = ''

	def __iter__(self):
		return self

	def next_read(self, reads):
		ID = reads.readline().strip().split(' ')[0]
		seq = reads.readline().strip()
		reads.readline()
		qual = reads.readline().strip()

		return ID, seq, qual

	def next(self):
		ID1, seq1, qual1 = self.next_read(self.reads1)
		while ID1 == self.last_ID1:
			ID1, seq1, qual1 = self.next_read(self.reads1)

		ID2, seq2, qual2 = self.next_read(self.reads2)
		while ID2 == self.last_ID2:
			ID2, seq2, qual2 = self.next_read(self.reads2)

		self.last_ID1, self.last_ID2 = ID1, ID2

		assert ID1 == ID2

		if ID1 == '':
			raise StopIteration

		#seq1 = revcomp(seq1)
		#qual1 = qual1[::-1]

		return ID1, seq1, seq2, qual1, qual2
