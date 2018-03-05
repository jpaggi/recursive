import pysam
from load_genome import *
import sys
MIN_INTRON_SIZE = 1000 

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
		chrom, start, end, strand = line.strip().split('\t')
		if strand == '+':
			five, three = int(start), int(end)
		else:
			three, five = int(start), int(end)
		chrom = chrom[3:]
		key = (chrom, strand, three)
		out[key] = out[key] + [five] if key in out else [five]
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
				three, five = blocks[i][1] - 1, blocks[i+1][0]

			# Some insertion events align to 5'ss and cause false positives
			if abs(five - three) < MIN_INTRON_SIZE: continue

			chrom = samfile.getrname(read.reference_id)
			strand_str = '+' if strand else '-'

			if (chrom, strand_str, three) in ss:
				intron_five = ss[(chrom, strand_str, three)]
				if strand and five > min(intron_five) and five - 1 not in intron_five and seq.query(chrom, strand, five):
					print '\t'.join(map(str, [chrom, five, three, read.query_name, blocks, strand_str]))
				elif not strand and five < max(intron_five) and five not in intron_five and seq.query(chrom, strand, five):
					print '\t'.join(map(str, [chrom, three, five+1, read.query_name, blocks, strand_str]))
