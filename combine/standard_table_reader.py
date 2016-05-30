CSV = lambda x: ','.join(map(str, list([x] if type(x) == int else x)))
#COORDS', 'RS', 'GRAV', 'SJR', 'SAW_SCORE', 'MCMC_PROB', 'MANUAL', 'SPANNING',
#	             'RECURSIVE_INDEX' 'MOTIF', 'MOTIF_SCORE', 'DOWN', 'BODY', 'FIVE_SCORE', 'PUT_FIVE_SCORE',
#	             'PUT_THREE', 'ANNO_THREE'
class Entry:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom, self.start, self.end, self.strand = self._read_coords(a[0])
		self.rs = int(a[1])
		self.grav = int(a[2])
		self.junc = map(int, a[3].split(','))
		self.saw_score = float(a[4])
		self.mcmc_prob = float(a[5])
		self.manual = int(a[6])
		self.intron_sjr = int(a[7])
		self.recursive_index = float(a[8])
		self.motif = a[9]
		self.motif_score = float(a[10])
		self.down_counts = map(int, a[11].split(','))
		self.body_counts = map(int, a[12].split(','))
		self.five_score  = float(a[13])
		self.put_five_score = float(a[14])
		self.putative_three = int(a[15])
		self.anno_three     = int(a[16])

	def _read_coords(self, coords):
		chrom, index, strand = coords.split(':')
		start, end = index.split('-')
		if strand == '+':
			start, end = set(map(int, start.split(','))), int(end)
		else:
			start, end = int(start), set(map(int, end.split(',')))
		return chrom, start, end, strand

	def _get_motif(self):
		motif = seq[self.chrom][self.rs - 20:self.rs + 20]
		if self.strand == '-': motif = revcomp(motif)
		return motif[:28]

	def set_intron_sjr(self, count):
		self.intron_sjr = count


	def good_motif(self, thresh):
		return self.motif_score > thresh

	def five(self):
		return self.start if self.strand =='+' else self.end

	def three(self):
		return self.end if self.strand =='+' else self.start

	# def add_sjr(self, counts):
	# 	self.junc = [i+j for i, j in zip(self.junc, counts)]

	# def add_sawtooth(self):
	# 	self.sawtooth = 1

	# def add_graveley(self):
	# 	self.grav = 1\

	def recursive_index(self):
		return sum(self.junc) / float(self.intron_sjr)

	def add_manual(self):
		self.manual = 1

	def decreasing_body_counts(self):
		diffs = [i - j for i, j in  zip(self.body_counts[1:], self.body_counts[:-1])]
		return len(filter(lambda x: x < 0, diffs)) > 2

	def _coords_str(self):
		return "{}:{}-{}:{}".format(self.chrom, CSV(self.start), CSV(self.end), self.strand)

	def five_min(self):
		return sum(self.junc[:2])

	def ten_min(self):
		return sum(self.junc[2:4])

	def twenty_min(self):
		return sum(self.junc[4:6])

	def total_min(self):
		return sum(self.junc[6:])

	def __str__(self):
		return '\t'.join(map(str, [
			self._coords_str(),
			self.rs,
			self.grav, 
			CSV(self.junc),
			self.saw_score,
			self.mcmc_prob,
			self.manual,
			self.intron_sjr,
			self.recursive_index,
			self.motif,
			self.motif_score,
			CSV(self.down_counts),
			CSV(self.body_counts),
			self.five_score,
			self.put_five_score,
			self.putative_three,
			self.anno_three]))

if __name__ == '__main__':
	import sys
	import matplotlib.pyplot as plt

	manual = []
	other = []

	for line in open(sys.argv[1]):
		if line[:6] == 'COORDS': continue
		entry = Entry(line)
		# if not entry.good_motif(.83): continue
		# if entry.manual:
		# 	print entry

		statistic = 2 * entry.total_min() + entry.twenty_min() - 2 * entry.five_min()  -  entry.ten_min()

		statistic = entry.motif_score

		if entry.manual:
			manual += [statistic]
		else:
			other += [statistic]

	print len(manual), len(filter(lambda x: x > .9, manual))

	print len(other), len(filter(lambda x: x > .9, other))

	plt.hist(manual, bins = 50)
	plt.axvline(0)
	plt.show()
	plt.hist(other, bins = 50)
	plt.axvline(0)
	plt.show()

