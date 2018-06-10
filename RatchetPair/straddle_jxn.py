import pysam
import sys
from load_genome import *
from seq import Seq

class Junctions:
	# Handles lookups of splice sites in a range of genomic coordinates.
	def __init__(self, file):
		# ss: strand -> chrom -> (five, three)
		self.ss = self._read_ss(file)
		self._sort()

	def _read_ss(self, file):
		ss = {'+': {}, '-': {}}
		with open(file) as sites:
			for line in sites:
				chrom, start, end, _, _, strand = line.strip().split('\t')
				if strand == '+':
					five, three = int(start), int(end)
				else:
					three, five = int(start), int(end)

				if chrom in ss[strand]:
					ss[strand][chrom] += [(five, three)]
				else:
					ss[strand][chrom] = [(five, three)]
		return ss

	def _sort(self):
		# Sort all strands, chromosomes by 5'ss
		for strand in ['+', '-']:
			for chrom in self.ss[strand]:
				self.ss[strand][chrom] = sorted(self.ss[strand][chrom], key=lambda x: x[0])

	def search(self, chrom, strand, start, end):
		if chrom not in self.ss['+'] or chrom not in self.ss['-']: return []
		up = self._linear_search(chrom, strand, start)
		down = self._linear_search(chrom, strand, end, up)
		return self.ss[strand][chrom][min(up, down):max(up, down)]

	def _linear_search(self, chrom, strand, pos, start = 0):
		for i in xrange(start, len(self.ss[strand][chrom])):
			if self.ss[strand][chrom][i][0] > pos:
				return i
		return len(self.ss[strand][chrom])

MAX_ALLOWED_INSERT = 1000
MAX_FIVE_INSERT = 400
OVERHANG = 10           # Don't extract reads short of an annotated splice site by less than this.
ALLOW_ANNO_AGGT = True  # Allow alignment to anno splice sites with seq AGGT and inside bigger intron.

samfile = pysam.AlignmentFile(sys.argv[1])
jxns    = Junctions(sys.argv[2])
seq = Seq(load_genome(open(sys.argv[3], 'r')))
sample = sys.argv[4]

for read in samfile.fetch():
	if not read.is_read2: continue
	if read.reference_id != read.next_reference_id: continue

	strand = (read.is_read1 == read.is_reverse)
	strand_str = '+' if strand else '-'
	chrom = samfile.getrname(read.reference_id)
	if chrom[:3] == 'chr': chrom = chrom[3:]

	# The following two blocks are identical and handle + and - strands.
	if strand and read.pnext - read.pos > MAX_ALLOWED_INSERT:
		inner_left = read.get_blocks()[-1][1] # End   of mate #1
		inner_right = read.pnext              # Start of mate #2
		
		# [(five, three), ...] for all splices between read pair
		splices = jxns.search(chrom, strand_str, inner_left, inner_right)

		# check if consistent with any anno splice
		consistent = False
		for five, three in splices:
			if ALLOW_ANNO_AGGT and seq.query(chrom, strand, three): continue
			cassette = five + OVERHANG >= inner_left and three - OVERHANG <= inner_right
			if cassette and five - inner_left + inner_right - three < MAX_ALLOWED_INSERT:
				consistent = True
		if consistent: continue

		fives = {}
		for five, three in splices:
			if five not in fives: fives[five] = three
			fives[five] = max(three, fives[five])

		# Only print if downstream end lies within an annotated intron
		# and 5'ss close to inner_left.
		for five, three in fives.items():
			if five - read.pos < MAX_FIVE_INSERT and inner_right + OVERHANG < three:
				print '\t'.join(map(str, [chrom, inner_left, inner_right, five, sample, strand_str]))

	elif not strand and MAX_ALLOWED_INSERT < read.pos - read.pnext:
		inner_left  = read.pnext
		inner_right = read.pos
		splices = jxns.search(chrom, strand_str, inner_left, inner_right)

		consistent = False
		for five, three in splices:
			if ALLOW_ANNO_AGGT and seq.query(chrom, strand, three): continue
			cassette = three + OVERHANG >= inner_left and five - OVERHANG <= inner_right
			if cassette and three - inner_left + inner_right - five < MAX_ALLOWED_INSERT:
				consistent = True
		if consistent: continue

		fives = {}
		for five, three in splices:
			if five not in fives: fives[five] = three
			fives[five] = min(three, fives[five])			

		for five, three in fives.items():
			if inner_right - five < MAX_FIVE_INSERT and inner_left > three + OVERHANG:
				print '\t'.join(map(str, [chrom, inner_left, inner_right, five, sample, strand_str]))
