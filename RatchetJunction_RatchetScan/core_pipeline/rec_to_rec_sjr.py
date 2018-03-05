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
	out5, out3 = {}, {}
	for line in sites:
		chrom, start, end, strand = line.strip().split('\t')
		if strand == '+':
			five, three = int(start)+1, int(end)
		else:
			three, five = int(start), int(end)-1
		chrom = chrom[3:]
		out3[(chrom, strand, three)] = 1
		out5[(chrom, strand, five)] = 1
		
	return out5, out3

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
ss5, ss3 = get_jxns(sys.argv[2])
seq = Seq(load_genome(open(sys.argv[3], 'r')))

for read in samfile.fetch():
	blocks = merge_blocks(read.get_blocks())
	strand = (read.is_read1 == read.is_reverse)
	if len(blocks) > 1:
		for i in xrange(len(blocks) - 1):
			if strand:
				five, three = blocks[i][1], blocks[i+1][0]
			else:
				three, five = blocks[i][1] - 1, blocks[i+1][0] - 1

			# Some insertion events align to 5'ss and cause false positives
			if abs(five - three) < MIN_INTRON_SIZE: continue

			chrom = samfile.getrname(read.reference_id)
			strand_str = '+' if strand else '-'

			if (chrom, strand_str, three) in ss3: continue
			if (chrom, strand_str, five) in ss5: continue

			if not strand: five, three = five + 1, three + 1
			if not seq.query(chrom, strand, five): continue
			if not seq.query(chrom, strand, three): continue
			
			if strand:
				print '\t'.join(map(str, [chrom, five+1, three, read.query_name, blocks, strand_str]))
			else:
				print '\t'.join(map(str, [chrom, three, five+1, read.query_name, blocks, strand_str]))
