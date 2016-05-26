
# Input will be table of introns with recursive sites

# set 5' ss and 3'ss

# upstream and downstream exons

import sys
import pysam

def merge_blocks(blocks):
	"""
	Removes short gaps in read that are likely indels
	"""
	i = 1
	while i < len(blocks) - 1:
		if blocks[i][0] - blocks[i-1][1] < 5:
			blocks[i] = (blocks[i][0], blocks[i+1][1])
			blocks.remove(blocks[i+1])
		else:
			i += 1
	return blocks

class Exon:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom = a[0]
		self.start = int(a[3])
		self.end   = int(a[4])
		self.strand= a[6]

def GTF_EXON(name, start, end, pos, mRNA):
	atts = "ID={}.mRNA.{}.exon.{};Parent={}.mRNA.{}".format(name, mRNA, pos, name, mRNA)
	return '\t'.join(map(str, [name, 'RI', 'exon', start, end, '.', '+', '.', atts]))

def GTF_MRNA(name, start, end, pos):
	atts = "ID={}.mRNA.{};Parent={}".format(name, pos, name)
	return '\t'.join(map(str, [name, 'RI', 'mRNA', start, end, '.', '+', '.', atts]))

def GTF_GENE(name, start, end):
	atts = "ID={}".format(name)
	return '\t'.join(map(str, [name, 'RI', 'gene', start, end, '.', '+', '.', atts]))

class Intron:
	"""
	Internally uses one-based coords

	Exon coords point to first or last base in exon

	RS coords point to AG[G]T
	"""
	def __init__(self, line):
		a = line.strip().split()
		self.chrom  = a[0]
		self.strand = a[5]
		self.start  = int(a[1])
		self.end    = int(a[2]) + 2
		self.name   = a[3]
		self.rs     = int(a[4])
		if self.strand == '+': self.rs += 1

		self.exon_start = None
		self.exon_end   = None

	def maybe_add_exon(self, exon):
		if self.chrom != exon.chrom: return
		if self.strand != exon.strand: return
		if exon.start == self.end:
			if self.strand == '+':
				end = exon.end
			else:
				end = exon.end - 1
			self.exon_end = min(self.exon_end, end) if self.exon_end else end
		if exon.end  == self.start:
			if self.strand == '+':
				start = exon.start + 2
			else:
				start = exon.start
			self.exon_start = max(self.exon_start, start)


	def index(self, pos):
		if self.strand == '+':
			return pos - self.exon_start + 1000
		else:
			return self.exon_end - pos + 1000  

	def _get_name(self):
		return "{}:{}-{}:{}".format(self.chrom, self.start, self.end, self.strand)

	def smush(self):
		coords = sorted(map(self.index, [self.exon_start, self.start, self.rs, self.end, self.exon_end]))
		up_start, up_end, rs, down_start, down_end = coords
		down_end   += rs - down_start
		down_start += rs - down_start
		return up_start, up_end, rs, down_start, down_end

	def get_annotations(self):
		up_start, up_end, rs, down_start, down_end = self.smush()
		# make gene
		gene_name = "{}".format(self._get_name())
		lines = [GTF_GENE(gene_name, up_start, down_end)]
		# make mRNAs
		lines += [GTF_MRNA(gene_name, up_start, down_end, 'in')]
		lines += [GTF_MRNA(gene_name, up_start, down_end, 'ex')]
		# make inclusion exon
		lines += [GTF_EXON(gene_name, up_start, down_end, 1, 'in')]
		# make exclusion exons
		lines += [GTF_EXON(gene_name, up_start, up_end, 1, 'ex')]
		lines += [GTF_EXON(gene_name, down_start, down_end, 2, 'ex')]

		return '\n'.join(lines)+'\n'

	def smush_down(self):
		coords = sorted(map(self.index, [self.exon_start, self.start, self.rs, self.end, self.exon_end]))
		up_start, up_end, rs, down_start, down_end = coords
		rs         += up_end - rs
		down_start += up_end - rs
		down_end   += up_end - rs
		return up_start, up_end, rs, down_start, down_end

	def get_annotations_down(self):
		up_start, up_end, rs, down_start, down_end = self.smush_down()
		# make gene
		gene_name = "{}.down".format(self._get_name())
		lines = [GTF_GENE(gene_name, up_start, down_end)]
		# make mRNAs
		lines += [GTF_MRNA(gene_name, up_start, down_end, 'in')]
		lines += [GTF_MRNA(gene_name, up_start, down_end, 'ex')]
		# make inclusion exon
		lines += [GTF_EXON(gene_name, up_start, down_end, 1, 'in')]
		# make exclusion exons
		lines += [GTF_EXON(gene_name, up_start, up_end, 1, 'ex')]
		lines += [GTF_EXON(gene_name, down_start, down_end, 2, 'ex')]

		return '\n'.join(lines)+'\n'


	def get_header(self):
		name = "{}".format(self._get_name())
		length = self.smush()[-1]

		return '\t'.join(['@SQ', 'SN:'+name, 'LN:'+str(length)])


	def get_header_down(self):
		name = "{}.down".format(self._get_name())
		length = self.smush()[-1]

		return '\t'.join(['@SQ', 'SN:'+name, 'LN:'+str(length)])


	def transfer_block(self, blocks):
		temp = []
		for block in blocks:
			temp += [(self.index(block[0])+1, self.index(block[1]))]
		if self.strand == '+':
			return temp
		out = []
		for block in temp[::-1]:
			out += [(block[1], block[0]-2)]
		return out


	def assign_reads(self, reads, out):
		coords = sorted(map(self.index, [self.exon_start, self.start, self.rs, self.end, self.exon_end]))
		up_start, up_end, rs, down_start, down_end = coords
		gene_name = "{}".format(self._get_name())
		for read in reads.fetch(self.chrom, self.exon_start, self.exon_end):
			blocks = self.transfer_block(merge_blocks(read.get_blocks()))
			if len(blocks) == 2:
				five, three = blocks[0][1], blocks[1][0]
				if five == up_end:
					if three == down_start:
						len2 = blocks[1][1] - blocks[1][0] + 1
						new_blocks = [blocks[0], (rs, rs + len2)]
					elif three == rs:
						new_blocks = blocks
					else:
						continue
				else:
					continue
			elif len(blocks) == 1:
				start, end = blocks[0]
				# check if read in upstream exon or rs segment
				if start < rs and end > up_start:
					new_blocks = blocks

				# check if read in downstream exon
				elif down_start < start < down_end:
					start -= down_start - rs
					end   -= down_start - rs
					new_blocks = [(start, end)]
			else:
				continue
			entry = SamEntry(read, self._get_name(), new_blocks)
			out.write(str(entry) + '\n')

	def assign_reads_down(self, reads, out):
		coords = sorted(map(self.index, [self.exon_start, self.start, self.rs, self.end, self.exon_end]))
		up_start, up_end, rs, down_start, down_end = coords
		gene_name = "{}.down".format(self._get_name())

		s_three = self.smush_down()[2]

		for read in reads.fetch(self.chrom, self.exon_start, self.exon_end):
			blocks = self.transfer_block(merge_blocks(read.get_blocks()))
			if len(blocks) == 2:
				five, three = blocks[0][1], blocks[1][0]
				if five == up_end:
					if three == down_start:
						len2 = blocks[1][1] - blocks[1][0] + 1
						new_blocks = [blocks[0], (s_three, s_three + len2)]
					elif three == rs:
						# convert into ungapped read !!!
						new_blocks = [(five, five + 50)]
					else:
						continue
				else:
					continue
			elif len(blocks) == 1:
				start, end = blocks[0]
				# check if read in upstream exon
				if up_start < end < up_end:
					new_blocks = blocks

				# check if read in rs segment or downstream exon
				elif start < down_end and end > rs:
					start -= rs - up_end
					end   -= rs - up_end
					new_blocks = [(start, end)]
			else:
				continue
			entry = SamEntry(read, self._get_name(), new_blocks)
			out.write(str(entry) + '\n')

	def __str__(self):
		return '\t'.join(map(str, [
			self.chrom,
			self.exon_start,
			self.start,
			self.rs,
			self.end,
			self.exon_end]))


class SamEntry:
	def __init__(self, read, name, blocks):
		self.blocks = blocks
		self.ID      = read.query_name
		self.flag    = 0
		self.chrom   = name
		self.pos     = max(0, blocks[0][0])
		self.mapq    =  255

		first_len = blocks[0][1] - self.pos
		if len(blocks) == 1:
			self.cigar = "{}M".format(51)
		else:
			assert len(blocks) == 2
			second_length = 51 - first_len
			intron_length = blocks[1][0] - blocks[0][1] - 1
			self.cigar = "{}M{}N{}M".format(first_len, intron_length, second_length)
		self.rnext   = '*'
		self.pnext   = 0
		self.tlen    =  0
		self.seq     = 'A' * 51
		self.quality = 'A' * 51

	def __str__(self):
		return '\t'.join(map(str, [
			self.ID,
			self.flag,
			self.chrom, 
			self.pos + 1, 
			self.mapq,
			self.cigar,
			self.rnext,
			self.pnext,
			self.tlen,
			self.seq,
			self.quality]))

# read in all introns
introns = [Intron(line) for line in open(sys.argv[1]) if len(line.split('\t')[4].split(',')) == 1]

# read in all exons
exons   = [Exon(line) for line in open(sys.argv[2])]

annotations = open('../data/test_annotations.gff', 'w')

# join exons to introns
for intron in introns:
	for exon in exons:
		intron.maybe_add_exon(exon)
	print intron

for intron in introns:
	#annotations.write(intron.get_annotations())
	annotations.write(intron.get_annotations_down())

samfile = pysam.AlignmentFile('../data/10_minute_aligning.bam')

out = open('../data/test_grav.sam', 'w')
out.write('\t'.join(['@HD', 'VN:1.0', 'SO:unsorted']) + '\n')
for intron in introns:
	#out.write(intron.get_header() + '\n')
	out.write(intron.get_header_down() + '\n')
out.write('\t'.join(['@PG', 'ID:joe', 'PN:joe', 'VN:1.0']) + '\n')
for intron in introns:
	#intron.assign_reads(samfile, out)
	intron.assign_reads_down(samfile, out)


# for each intron write annotations

# for each BAM file
# .... for each intron write to new samfile
