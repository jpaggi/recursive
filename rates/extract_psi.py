import sys

"""
Creates a table of all miso mean psi values + 
identiy of intron and recursive segment being considered
"""

SAMPLES = ['5_repA', '5_repB', '5_repC', '10_repA', '10_repB', '10_repC',
	       '20_repA', '20_repB', '20_repC', 'total_repA', 'total_repB']
INDEX = lambda x: SAMPLES.index(x)
MERGE = lambda x, y: [max(i, j) for i, j in zip(x, y)]

BASE = sys.argv[1]
NAME = 'miso_summary/summary/miso.miso_summary'

CSV = lambda x: ','.join(map(str, x))


# example miso summary line
#X:2536126-2563946:+.rs1	0.03	0.02	0.03	'X:2536126-2563946:+.rs1.mRNA.in.exon.1','X:2536126-2563946:+.rs1.mRNA.ex.exon.1_X:2536126-2563946:+.rs1.mRNA.ex.exon.2'
#	(0,0):149,(0,1):77,(1,0):725,(1,1):144	0:728,1:218	X:2536126-2563946:+.rs1	+	1000,1000	17447,17447

class Intron:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom  = a[0]
		self.strand = a[5]
		self.start  = int(a[1])
		self.end    = int(a[2]) + 2
		self.name   = a[3]
		self.rs     = sorted(map(lambda x: int(x) + 1 if self.strand == '+' else int(x), a[4].split(',')))

		self.miso = []

	def key(self):
		return "{}:{}-{}:{}".format(self.chrom, self.start, self.end, self.strand)

	def add_miso(self, miso):
		self.miso += [miso]

	def get_rs_coords(self, i):
		bounds = [self.start] + self.rs + [self.end]
		if self.strand == '-': bounds = bounds[::-1]
		coords = bounds[i:i+2]
		if self.strand == '-': coords = coords[::-1]
		return coords

	def __str__(self):
		out = []
		for miso in self.miso:
			rs_start, rs_end = self.get_rs_coords(miso.num())
			out += [map(str, [
				self.chrom,
				self.start,
				self.end,
				self.strand,
				rs_start,
				rs_end,
				miso.num(),
				len(self.rs) + 1,
				CSV(miso.mean_psi)])]
		out = sorted(out, key = lambda x: x[6])
		return '\n'.join(map('\t'.join, out))

class MISO:
	def __init__(self, line, sample):
		a = line.strip().split()

		self.mean_psi = [-1] * len(SAMPLES)
		self.low_psi  = [-1] * len(SAMPLES)
		self.high_psi = [-1] * len(SAMPLES)
		self.counts   = ['*'] * len(SAMPLES)

		self.gene = a[0]
		self.mean_psi[INDEX(sample)]  = float(a[1])
		self.low_psi[INDEX(sample)]   = float(a[2])
		self.high_psi[INDEX(sample)]  = float(a[3])
		self.counts[INDEX(sample)]    = a[5]

	def key(self):
		return self.gene

	def merge(self, other):
		self.mean_psi  = MERGE(self.mean_psi, other.mean_psi)
		self.low_psi   = MERGE(self.low_psi, other.low_psi)
		self.high_psi  = MERGE(self.high_psi, other.high_psi)
		self.counts    = [i if j == '*' else j for i, j in zip(self.counts, other.counts)]

	def num(self):
		return int(self.gene.split('rs')[1])

	def __str__(self):
		return '\t'.join(map(str, [
			self.gene,
			CSV(self.mean_psi),
			CSV(self.low_psi),
			CSV(self.high_psi),
			';'.join(self.counts)]))


entries = {}
for sample in SAMPLES:
	f = open("{}/{}/{}".format(BASE, sample, NAME), 'r')
	f.readline() # pass the header...
	for line in f:
		miso = MISO(line, sample)
		if miso.key() in entries:
			entries[miso.key()].merge(miso)
		else:
			entries[miso.key()] = miso

print "#" + ','.join(SAMPLES)
print '\t'.join(['chrom', 'start', 'end', 'rs_start', 'rs_end', 'rs_num', 'total_rs', 'mean_psi'])


for line in open(sys.argv[2], 'r'):
	intron = Intron(line)

	for entry in entries:
		if intron.key() in entry:
			intron.add_miso(entries[entry])

	if intron.miso:
		print intron




