
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
		self.rs     = map(lambda x: int(x) + 1, a[4].split(','))

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

	def _get_name(self, i):
		return "{}:{}-{}:{}.rs{}".format(self.chrom, self.start, self.end, self.strand, i)

	def transfer_coords(self):
		return sorted(map(self.index, [self.exon_start, self.start] + self.rs + [self.end, self.exon_end]))

	def smush(self, i):
		coords = self.transfer_coords()
		up_start, up_end = coords[:2]
		down_start, down_end = coords[-2:]
		rs_start, rs_end = coords[1+i:3+i]

		rs_end    -= rs_start - up_end

		down_start -= down_start - rs_end
		down_end   -= down_start - rs_end
		return up_start, up_end, down_start, down_end

	def get_annotations(self, i):
		up_start, up_end, down_start, down_end = self.smush(i)
		# make gene
		gene_name = self._get_name(i)
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


	def get_header(self, i):
		name = self._get_name(i)
		length = self.smush(i)[-1] + 1000

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


	def assign_reads(self, i, reads, out):
		coords = self.transfer_coords()
		up_start, up_end = coords[:2]
		down_start, down_end = coords[-2:]

		sites = coords[1:-1]

		rs_start, rs_end = sites[i:i+2]
		exc = sites[i+1:]

		s_three = self.smush(i)[2]

		gene_name = self._get_name(i)
		for read in reads.fetch(self.chrom, self.exon_start, self.exon_end):
			blocks = self.transfer_block(merge_blocks(read.get_blocks()))
			if len(blocks) == 2:
				five, three = blocks[0][1], blocks[1][0]
				if five == up_end:
					if three in exc:
						len2 = blocks[1][1] - blocks[1][0] + 1
						new_blocks = [blocks[0], (s_three, s_three + len2)]
					elif three == rs_start:
						new_blocks = [(five, five + 50)]
					else:
						continue
				else:
					continue
			elif len(blocks) == 1:
				start, end = blocks[0]
				# check if read entirely in upstream exon
				if up_start < end < up_end:
					new_blocks = blocks
				# check if read entirely in downstream exon
				elif down_start < start < down_end:
					start -= down_start - s_three
					end   -= down_start - s_three
					new_blocks = [(start, end)]
				# check if overlaps recursive segment
				elif start < rs_end and end > rs_start:
					start -= rs_start - up_end
					end   -= rs_start - up_end
					new_blocks = [(start, end)]
				else:
					continue
			else:
				continue
			entry = SamEntry(read, self._get_name(i), new_blocks)
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
introns = [Intron(line) for i, line in enumerate(open(sys.argv[1])) if i < 20]

# read in all exons
exons   = [Exon(line) for line in open(sys.argv[2])]

annotations = open('../data/test_annotations.gff', 'w')

# join exons to introns
for intron in introns:
	for exon in exons:
		intron.maybe_add_exon(exon)
	print intron
	print intron.rs

for intron in introns:
	for i in range(len(intron.rs)+1):
		annotations.write(intron.get_annotations(i))

samfile = pysam.AlignmentFile('../data/10_minute_aligning.bam')

out = open('../data/test_grav.sam', 'w')

out.write('\t'.join(['@HD', 'VN:1.0', 'SO:unsorted']) + '\n')
for intron in introns:
	for i in range(len(intron.rs)+1):
		out.write(intron.get_header(i) + '\n')
out.write('\t'.join(['@PG', 'ID:joe', 'PN:joe', 'VN:1.0']) + '\n')

for intron in introns:
	for i in range(len(intron.rs)+1):
		intron.assign_reads(i, samfile, out)


# for each intron write annotations

# for each BAM file
# .... for each intron write to new samfile
