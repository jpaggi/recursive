CSV = lambda x: ','.join(map(str, list([x] if type(x) == int else x)))

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
		self.log_or      = float(a[11])
		self.down_counts = map(int, a[12].split(','))
		self.body_counts = map(int, a[13].split(','))
		self.five_score  = float(a[14])
		self.put_five_score = float(a[15])
		self.putative_three = int(a[16])
		self.putative_five  = int(a[17])
		self.anno_five     = int(a[18])

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

	def set_putative(self, five, three, score):
		self.putative_five = five
		self.putative_three = three
		self.put_five_score = score

	def good_motif(self, thresh):
		return self.motif_score > thresh

	def five(self):
		return self.start if self.strand =='+' else self.end

	def three(self):
		return self.end if self.strand =='+' else self.start

	def recursive_index(self):
		return sum(self.junc) / float(self.intron_sjr)

	def add_manual(self):
		self.manual = 1

	def decreasing_body_counts(self):
		diffs = [j-i for i, j in  zip(self.body_counts[1:], self.body_counts[:-1])]
		return len(filter(lambda x: x >= 0, diffs)) > 1

	def decreasing_junc_counts(self):
		counts = [self.five_min(), self.ten_min(), self.twenty_min(), self.total_min()]
		diffs = [j-i for i, j in  zip(counts[1:], counts[:-1])]
		return len(filter(lambda x: x >= 0, diffs)) > 1

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
			self.log_or,
			CSV(self.down_counts),
			CSV(self.body_counts),
			self.five_score,
			self.put_five_score,
			self.putative_three,
			self.putative_five,
			self.anno_five]))

if __name__ == '__main__':
	import sys
	count = 0
	for line in open(sys.argv[1]):
		if line[:6] == 'COORDS': continue
		entry = Entry(line)
		if entry.motif[23:27] != 'AGGT':
			print entry.junc
			count += 1

		sjr = entry.recursive_index > .05 and sum(entry.junc) and entry.log_or > 0
		saw = entry.saw_score > .06
		exon = sum(entry.down_counts) * 5 > sum(entry.junc) and not entry.decreasing_body_counts()

		if not sjr and not saw: continue
		if exon: continue
		if entry.grav: continue

		if entry.manual != -1: continue

		print entry
	print count
