
def get_jxns(anno_file):
	"""
	Return a dictionary of splicing events present in the anno_file

	The anno file should have exons named well 'exon' and can be in
	either gff or gtf format

	The dictionary is (chrom, strand, five) -> [<associated 3'ss>, ...]
	"""
	
	sites =  open(anno_file, 'r')

	out = {}
	for line in sites:
		chrom, start, end, strand = line.strip().split('\t')[:7]
		if strand == '+':
			five, three = int(start), int(end)
		else:
			three, five = int(start), int(end)

		key = (chrom, strand, five)

		if key in out:
			out[key].append(three)
		else:
			out[key] = [three]
	return out


class Junctions:
	def __init__(self, file):
		# strand -> chrom -> (five, three)
		self.ss = {'+': {}, '-': {}}
		sites = open(file, 'r')
		for line in sites:
			chrom, start, end, strand = line.strip().split('\t')[:7]
			if strand == '+':
				five, three = int(start), int(end)
			else:
				three, five = int(start), int(end)

			if chrom in self.ss[strand]:
				self.ss[strand][chrom] += [(five, three)]
			else:
				self.ss[strand][chrom] = [(five, three)]
		self.sort()
		
	def sort(self):
		for strand in ['+', '-']:
			for chrom in self.ss[strand]:
				self.ss[strand][chrom] = sorted(self.ss[strand][chrom], key=lambda x: x[0])


	def _linear_search(self, chrom, strand, pos, start = 0):
		for i in xrange(start, len(self.ss[strand][chrom])):
			if self.ss[strand][chrom][i][0] > pos:
				return i
		return len(self.ss[strand][chrom])


	def search(self, chrom, strand, start, end):
		up = self._linear_search(chrom, strand, start)
		down = self._linear_search(chrom, strand, end, up)
		return [] if up == down else self.ss[strand][chrom][min(up, down):max(up, down)]
