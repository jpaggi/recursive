import sys
import pysam

def merge_blocks(blocks):
	"""
	Removes short gaps in read that are likely indels
	"""
	i = 1
	while i < len(blocks) - 1:
		if blocks[i][0] - blocks[i-1][1] < 10:
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

		down_end   -= down_start - rs_end
		down_start -= down_start - rs_end
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


	def assign_reads(self, i, reads):
		coords = self.transfer_coords()
		up_start, up_end = coords[:2]
		down_start, down_end = coords[-2:]
		gene_name = self._get_name(i)

		sites = coords[1:-1]

		rs_start, rs_end = sites[i:i+2]
		exc = sites[i+1:]

		s_three = self.smush(i)[2]

		entries = {}
		for read in reads.fetch(self.chrom, self.exon_start, self.exon_end):
			if read.is_unmapped: continue
			blocks = self.transfer_block(merge_blocks(read.get_blocks()))
			conditional = False
			new_blocks = None
			for k in range(len(blocks) -1):
				five, three = blocks[k][1], blocks[k+1][0]
				if five == up_end:
					up_length = sum(blocks[j][1] - blocks[j][0] + 1 for j in range(k+1))
					if three in exc:
						new_blocks = [(five - up_length, up_end), (s_three, s_three + 50 - up_length)]
					elif three == rs_start:
						new_blocks = [(five - up_length, five - up_length + 50)]

			if not new_blocks:
				start, end = blocks[0][0], blocks[-1][1]
				# check if read entirely in upstream exon
				if up_start < end < up_end:
					new_blocks = [(end - 50, end)]
					conditional = exc[0] < self.index(read.pos) < down_start
				# check if read entirely in downstream exon
				elif down_start < start < down_end:
					start -= down_start - s_three
					new_blocks = [(start, start+50)]
				# check if overlaps recursive segment

			if not new_blocks and len(blocks) == 1:
				if start < rs_end and end > rs_start:
					start -= rs_start - up_end
					new_blocks = [(start, start+50)]
				elif any(0 < start-p < 300 for p in exc) and self.index(read.pos) < up_end:
					base = [0 < start-p < 300 for p in exc][0]
					start -= base - s_three
					new_blocks = [(start, start+50)]
					conditional = True

			if new_blocks:
				if read.query_name in entries:
					entries[read.query_name].add(SamEntry(read, self._get_name(i), new_blocks))
				else:
					entries[read.query_name] = Pair(SamEntry(read, self._get_name(i), new_blocks), conditional, self)

		return filter(lambda x: x.is_valid(), entries.values())


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
		"""
		Initialize as being single end aligning to reverse strand
		"""
		self.blocks  = blocks
		self.ID      = read.query_name
		self.flag    = 1 + 8 + 64 # mulitple segs, mate unmapped, reversed, first
		self.chrom   = name
		self.pos     = max(0, blocks[0][0])
		self.mapq    =  255

		first_len = blocks[0][1] - self.pos
		if len(blocks) == 1:
			self.cigar = '51M'
			self.end = self.pos + 51
		else:
			assert len(blocks) == 2
			second_length = 51 - first_len
			intron_length = blocks[1][0] - blocks[0][1] - 1
			self.cigar = "{}M{}N{}M".format(first_len, intron_length, second_length)
			self.end = self.pos + 51 + intron_length

		self.rnext   = '='
		self.pnext   = self.pos
		self.tlen    =  0
		self.seq     = 'A' * 51
		self.quality = 'A' * 51

	def phantom_mate(self):
		return '\t'.join(map(str, [
			self.ID,
			133,
			self.chrom,
			self.pos,
			0,
			'*',
			'=',
			self.pos,
			0,
			self.seq,
			self.quality]))

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

class Pair:
	def __init__(self, sam_entry, conditional,  intron):
		self.read1 = sam_entry
		self.read2 = None
		self.conditional = conditional
		self.intron = intron

	def add(self, sam_entry):
		if sam_entry.pos > self.read1.pos:
			self.read2 = sam_entry
		else:
			self.read1, self.read2 = sam_entry, self.read1

		# if conditional check if good
		if self.conditional:
			# at the moment just getting the other read through is enough!
			self.conditional = False

		# set tlen
		self.read1.tlen = self.read2.end - self.read1.pos
		self.read2.tlen = - self.read1.tlen

		# set pnext
		self.read1.pnext = self.read2.pos
		self.read2.pnext = self.read1.pos

		# set rnext
		self.read1.rnext = self.read1.chrom
		self.read2.rnext = self.read2.chrom

		# set flag
		self.read1.flag = 1 + 2 + 64 + 32  # multiple, proper, first, reverse
		self.read2.flag = 1 + 2 + 128 + 16 # multiple, proper, first, mate reverse

	def is_valid(self):
		return not self.conditional

	def __str__(self):
		if self.read2:
			return str(self.read1) + '\n' + str(self.read2)
		else:
			return str(self.read1) + '\n' + self.read1.phantom_mate()

# read in all introns
introns = [Intron(line) for i, line in enumerate(open(sys.argv[1])) if i < 20]

# read in all exons
exons   = [Exon(line) for line in open(sys.argv[2])]

annotations = open("{}/annotations.gff".format(sys.argv[4]), 'w')

# join exons to introns
for intron in introns:
	for exon in exons:
		intron.maybe_add_exon(exon)
	print intron
	print intron.rs

for intron in introns:
	for i in range(len(intron.rs)+1):
		annotations.write(intron.get_annotations(i))



samfile = pysam.AlignmentFile(sys.argv[3])

out = open("{}/reads.sam".format(sys.argv[4]), 'w')

out.write('\t'.join(['@HD', 'VN:1.0', 'SO:unsorted']) + '\n')
for intron in introns:
	for i in range(len(intron.rs)+1):
		out.write(intron.get_header(i) + '\n')
out.write('\t'.join(['@PG', 'ID:joe', 'PN:joe', 'VN:1.0']) + '\n')

for intron in introns:
	for i in range(len(intron.rs)+1):
		out.write('\n'.join(map(str, intron.assign_reads(i, samfile))) + '\n')

