import numpy as np

class Repeats:
	def __init__(self, file):
		self.repeats = {}
		with open(file, 'r') as fp:
			for line in fp:
				if 'chr' not in line: continue
				chrom, start, end = line.strip().split()[4:7]
				start, end, chrom = int(start), int(end), chrom[3:]
				if chrom in self.repeats:
					self.repeats[chrom] += [(start, end)]
				else:
					self.repeats[chrom] = [(start, end)]
		for chrom in self.repeats:
			self.repeats[chrom] = sorted(self.repeats[chrom], key = lambda x: x[0])

	def _linear_search(self, chrom, pos):
		for i, (start, end) in enumerate(self.repeats[chrom]):
			if start >= pos:
				return i
		return len(self.repeats[chrom])

	def get_repeats(self, chrom, start, end):
		s, e = self._linear_search(chrom, start), self._linear_search(chrom, end)
		r = []
		for begin, end in self.repeats[chrom][s:e]:
			r += [(begin, end)]
		return r

	def mask(self, chrom, start, end, strand, expression):
		for b, s in self.get_repeats(chrom, start, end):
			if strand == '+':
				begin, stop = b - start, s - start
			else:
				begin, stop = end - s, end - b
			left = expression[max(0, begin - 1000): begin - 100]
			right = expression[stop + 100: min(len(expression), stop + 1000)]
			if left + right:
				median = int(np.median(np.asarray(left + right)))
			else:
				median = 0
			for pos in range(max(0, begin-100), min(stop+100, len(expression))):
				expression[pos] = median
		return expression
