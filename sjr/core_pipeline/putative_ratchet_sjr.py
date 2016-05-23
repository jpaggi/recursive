import pysam
from load_genome import *
import sys
MIN_INTRON_SIZE = 35

"""
Finds all reads in the given samfile that are broken and have a
5'ss aligning to a known 5'ss and a 3'ss aligning upstream of
a known 3'ss for the given intron. And have an AGGT at 3' end.

Constrains that abs(5' - 3') is greater than MIN_INTRON_SIZE.

Need to tune expression for strand based on geometry of library.
"""

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

def get_jxns(anno_file):
	sites =  open(anno_file, 'r')
	out = {}
	for line in sites:
		chrom, start, end, name, score, strand = line.strip().split('\t')[:6]
		if strand == '+':
			five, three, compare = int(start), int(end), max
		else:
			three, five, compare = int(start), int(end), min

		key = (chrom, strand, five)
		out[key] = compare(out[key], three) if key in out else three
	return out

class Seq:
	def __init__(self, seq):
		self.seq = seq
		self.AGGT = {}

	def query(self, chrom, strand, position):
		key = (chrom, strand, position)
		if not key in self.AGGT:
			self._add_to_AGGT(chrom, strand, position)
		return self.AGGT[key]

	def _add_to_AGGT(self, chrom, strand, position):
		key = (chrom, strand, position)
		motif = self.seq[chrom][position-2:position+2]
		if not strand: motif = revcomp(motif)
		self.AGGT[key] = (motif == 'AGGT')


samfile = pysam.AlignmentFile(sys.argv[1], 'rb')
ss = get_jxns(sys.argv[2])
seq = Seq(load_genome(open(sys.argv[3], 'r')))

for read in samfile.fetch():
	blocks = merge_blocks(read.get_blocks())
	strand = (read.is_read1 == read.is_reverse)
	if len(blocks) > 1:
		for i in xrange(len(blocks) - 1):
			if strand:
				five, three = blocks[i][1], blocks[i+1][0]
			else:
				three, five = blocks[i][1], blocks[i+1][0] - 1

			# Some insertion events align to 5'ss and cause false positives
			if abs(five - three) < MIN_INTRON_SIZE: continue

			chrom = samfile.getrname(read.reference_id)
			strand_str = '+' if strand else '-'

			if (chrom, strand_str, five) in ss:
				intron_three = ss[(chrom, strand_str, five)]
				if strand and three < intron_three and seq.query(chrom, strand, three):
					print '\t'.join(map(str, [chrom, five, three, read.query_name, blocks, strand_str]))
				elif not strand and three > intron_three and seq.query(chrom, strand, three):
					print '\t'.join(map(str,[chrom, three, five+1, read.query_name, blocks, strand_str]))
