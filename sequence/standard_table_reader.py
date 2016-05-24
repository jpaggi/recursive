CSV = lambda x: ','.join(map(str, list([x] if type(x) == int else x)))

class Entry:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom, self.start, self.end, self.strand = self._read_coords(a[0])
		self.rs = int(a[1])
		self.grav = int(a[2])
		self.junc = map(int, a[3].split(','))
		self.sawtooth = int(a[4])
		self.manual = int(a[5])
		self.motif_score = float(a[6])
		self.intron_sjr = int(a[7])

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
	# 	self.grav = 1

	def add_manual(self):
		self.manual = 1

	def _coords_str(self):
		return "{}:{}-{}:{}".format(self.chrom, CSV(self.start), CSV(self.end), self.strand)

	def __str__(self):
		return '\t'.join(map(str, [
			self._coords_str(),
			self.rs,
			self.grav,
			CSV(self.junc),
			self.sawtooth,
			self.manual,
			self.motif_score,
			self.intron_sjr]))

if __name__ == '__main__':
	import sys
	for line in open(sys.argv[1]):
		if line[:6] == 'COORDS': continue
		entry = Entry(line)
		if not entry.good_motif(.83): continue
		if entry.manual and (sum(entry.junc) or entry.grav):
			print entry

