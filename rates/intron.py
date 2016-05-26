
# Input will be table of introns with recursive sites

# set 5' ss and 3'ss

# upstream and downstream exons

import sys

class Exon:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom = a[0]
		self.start = int(a[3])
		self.end   = int(a[4])
		self.strand= a[6]

def GTF_EXON(name, start, end, pos, mRNA):
	atts = "ID={}.exon.{},Parent={}.mRNA.{}".format(name, pos, name, mRNA)
	return '\t'.join(map(str, [name, 'RI', 'exon', start, end, '.', '+', '.', atts]))

def GTF_MRNA(name, start, end, pos):
	atts = "ID={}.mRNA.{},Parent={}".format(name, pos, name)
	return '\t'.join(map(str, [name, 'RI', 'mRNA', start, end, '.', '+', '.', atts]))

def GTF_GENE(name, start, end):
	atts = "ID={}".format(name)
	return '\t'.join(map(str, [name, 'RI', 'gene', start, end, '.', '+', '.', atts]))

class Intron:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom  = a[0]
		self.start  = int(a[1])
		self.end    = int(a[2])
		self.name   = a[3]
		self.rs     = sorted(map(int, a[4].split(',')))
		self.strand = a[5]

		self.exon_start = None
		self.exon_end   = None

	def maybe_add_exon(self, exon):
		if self.chrom != exon.chrom: return
		if self.strand != exon.strand: return
		if exon.start == self.end + 2:
			self.exon_end = min(self.exon_end, exon.end) if self.exon_end else exon.end
		if exon.end  == self.start:
			self.exon_start = max(self.exon_start, exon.start)


	def index(self, pos):
		if self.strand == '+':
			return pos - self.exon_start + 1
		else:
			return self.exon_end - pos + 1

	def _get_name(self):
		return "{}:{}-{}:{}".format(self.chrom, self.start, self.end, self.strand)

	def length(self):
		return self.exon_end - self.exon_start

	def get_segment(self, i):
		sites = sorted([self.start] + self.rs + [self.end])
		if self.strand == '-': sites = sites[::-1]
		return sites[i:i+2]

	def get_exons(self):
		coords = map(self.index, [self.exon_start, self.start, self.exon_end, self.end])
		return sorted(coords)

	def get_transfer(self, i):
		rs_start, rs_end = map(self.index, self.get_segment(i))
		up_start, up_end, down_start, down_end = self.get_exons()

		gap1 = lambda x: x - rs_start + up_end
		gap2 = lambda x: x - down_start + rs_end

		regions = {(up_start, up_end): lambda x: x,
				   (rs_start, rs_end): gap1,
				   (down_start, down_end): lambda x: gap2(gap1(x))}
		def transfer(x):
			x = self.index(x)
			for region in regions:
				if region[0] <= x <= region[1]:
					return regions[region](x)
			return 0
		return transfer

	def get_annotations(self, i):
		transfer = self.get_transfer(i)
		coords = sorted(map(transfer, [self.exon_start, self.start] + self.get_segment(i) + [self.end, self.exon_end]))
		up_start, up_end, rs_start, rs_end, down_start, down_end = coords
		# make gene
		gene_name = "{}.rs{}".format(self._get_name(), i)
		lines = [GTF_GENE(gene_name, up_start, down_end)]
		# make mRNAs
		lines += [GTF_MRNA(gene_name, up_start, down_end, 'in')]
		lines += [GTF_MRNA(gene_name, up_start, down_end, 'ex')]
		# make inclusion exon
		lines += [GTF_EXON(gene_name, up_start, down_end, 1, 'in')]
		# make exclusion exons
		lines += [GTF_EXON(gene_name, up_start, up_end, 1, 'ex')]
		lines += [GTF_EXON(gene_name, down_start, down_end, 1, 'ex')]

		return '\n'.join(lines).strip()


	def assign_reads(self, i, reads):
		coords = sorted(map(self.index, [self.exon_start, self.start] + self.get_segment(i) + [self.end, self.exon_end]))
		up_start, up_end, rs_start, rs_end, down_start, down_end = coords

		gene_name = "{}.rs{}".format(self._get_name(), i)
		for read in reads.fetch(self.chrom, self.exon_start, self.exon_end):
			entry = SamEntry(read)

			# check if read in upstream exon

			# check if read in downstream exon

			# check if read is in recursive segment

			# check if read is junction from upstream to rs

			# check if read is junction into rs from upstream rs


class SamEntry:
	def __init(self, read, name):
		self.ID      = read.query_name
		self.flag    = read.FLAG
		self.chrom   = name
		self.pos     = None
		self.mapq    = read.mapq
		self.cigar   = None
		self.rnext   = name
		self.pnext   = None
		self.tlen    = None
		self.seq     = '.'
		self.quality = '.'

		self.blocks = read.get_blocks()

	def set_pos(self, pos):
		self.pos = pos

	def __str__(self):
		return '\t'.join(map(str, [
			self.ID,
			self.flag,
			self.chrom, 
			self.pos, 
			self.mapq,
			self.cigar,
			self.rnext,
			self.pnext,
			self.tlen,
			self.seq,
			self.quality]))

# read in all introns
introns = [Intron(line) for line in open(sys.argv[1])]

# read in all exons
exons   = [Exon(line) for line in open(sys.argv[2])]


# join exons to introns
for intron in introns:
	for exon in exons:
		intron.maybe_add_exon(exon)

for intron in introns:
	for i in range(len(intron.rs)+1):
		print intron.get_annotations(i)

# for each intron write annotations

# for each BAM file
# .... for each intron write to new samfile
