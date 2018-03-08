from anno_jxns import Junctions
import pysam
import sys
from load_genome import *

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

samfile = pysam.AlignmentFile(sys.argv[1])
jxns    = Junctions(sys.argv[2])
sample  = sys.argv[3]
seq = Seq(load_genome(open(sys.argv[4], 'r')))

for read in samfile.fetch():
	if not read.is_read2: continue
	if read.reference_id != read.next_reference_id: continue

	strand = (read.is_read1 == read.is_reverse)
	strand_str = '+' if strand else '-'
	chrom = samfile.getrname(read.reference_id)
	if chrom == 'Uextra': continue
	if chrom == 'dmel_mitochondrion_genome': continue

	if strand and 1000 < read.pnext - read.pos:
		inner_left = read.get_blocks()[-1][1]
		inner_right = read.pnext
		splices = jxns.search('chr' + chrom, strand_str, inner_left, inner_right)

		# check if consistent with any anno splice
		consistent = False
		for five, three in splices:
			if seq.query(chrom, strand, three): continue
			cassette = five + 10 >= inner_left and three - 10 <= inner_right
			if cassette and five - inner_left + inner_right - three < 600:
				consistent = True
		if consistent: continue

		fives = {}
		for five, three in splices:
			if five in fives:
				fives[five] = max(three, fives[five])
			else:
				fives[five] = three

		for five in fives:
			three = fives[five]
			if five - read.pos < 400 and inner_right + 10 < three:
				print '\t'.join(map(str, [chrom, inner_left, inner_right, sample, five + 1, strand_str]))

	elif not strand and 1000 < read.pos - read.pnext:
		inner_left  = read.pnext
		inner_right = read.pos
		splices = jxns.search('chr' + chrom, strand_str, inner_left, inner_right)

		consistent = False
		for five, three in splices:
			if seq.query(chrom, strand, three+1): continue
			cassette = three + 10 >= inner_left and five - 10 <= inner_right
			if cassette and three - inner_left + inner_right - five < 600:
				consistent = True
		if consistent: continue

		fives = {}
		for five, three in splices:
			if five in fives:
				fives[five] = min(three, fives[five])
			else:
				fives[five] = three

		for five in fives:
			three = fives[five]
			if inner_right - five < 400 and inner_left > three + 10:
				print '\t'.join(map(str, [chrom, inner_left, inner_right, sample, five, strand_str]))