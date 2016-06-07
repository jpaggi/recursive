
# Iterate through putative recursive sites...

# 	1) Iterate though sjr in 500 downstream
# 	2) Extract any sjr
# 	2) Fetch sequence 500 downstream
# 	3) search for all 5' ss motifs above a given threshold

import sys
import pysam
from standard_table_reader import Entry
from load_genome import *
from get_motifs import *
MAX_EXON = 500
MIN_MOTIF_SCORE = .5

def merge_blocks(blocks):
    """
    Removes short gaps in read that are likely indels
    """
    i = 1
    while i < len(blocks) - 1:
        if blocks[i][0] - blocks[i-1][1] < INDEL_LEN:
            blocks[i] = (blocks[i][0], blocks[i+1][1])
            blocks.remove(blocks[i+1])
        else:
            i += 1
    return blocks

class RS:
	def __init__(self, entry, seq):
		self.rs = entry
		self.sjr = {}  # position --> count???

		self.seq = seq

		self.three = None


	def _score_five(self):
		five = self.three[0]
		if not five: return 0
		if self.rs.strand == '+':
			begin = five
			end   = five + len(fp_pwm)
			motif = self.seq[self.rs.chrom][begin:end]
		else:
			begin = five - len(fp_pwm)
			end = five
			motif =  revcomp(self.seq[self.rs.chrom][begin:end])

		return (score_motif(fp_pwm, motif) - f_min) / (f_max - f_min)

	# def _set_seq(self, seq):
	# 	if self.rs.strand == '+':
	# 		begin = self.rs.rs
	# 		end   = self.rs.rs + len(fp_pwm) + MAX_EXON
	# 		return seq[self.rs.chrom][begin:end]
	# 	else:
	# 		begin = self.rs.rs - len(fp_pwm) - MAX_EXON
	# 		end = self.rs.rs
	# 		return revcomp(seq[self.rs.chrom][begin:end])

	# def _best_motif(self):
	# 	best = (-1, MIN_MOTIF_SCORE)
	# 	for i in range(1, len(self.seq) - len(fp_pwm)):  # +1 to not score the putative rs :)
	# 		score = (score_motif(fp_pwm, self.seq[i:i+len(fp_pwm)]) - f_min) / (f_max - f_min)
	# 		if score > best[1]:
	# 			best = (i, score)
	# 	if best[0] == -1:
	# 		return '.'
	# 	elif self.rs.strand == '+':
	# 		return self.rs.rs + best[0]
	# 	else:
	# 		return self.rs.rs - best[0]

	def add_sjr(self, position, three):
		if self._in_range(position):
			if (position, three) in self.sjr:
				self.sjr[(position, three)] += 1
			else:
				self.sjr[(position, three)] = 1

	def _assign_three(self):
		if self.sjr:
			self.three = [x for x in self.sjr if self.sjr[x] == max(self.sjr.values())][0]
		else:
			self.three = (0, 0)

	def _in_range(self, position):
		if self.rs.strand == '+':
			begin = self.rs.rs
			end   = self.rs.rs + MAX_EXON
		else:
			begin = self.rs.rs - MAX_EXON
			end = self.rs.rs

		return begin < position < end

	def get_coords(self):
		if self.rs.strand == '+':
			begin = self.rs.rs
			end   = self.rs.rs + MAX_EXON
		else:
			begin = self.rs.rs - MAX_EXON
			end = self.rs.rs

		return (self.rs.chrom, begin, end)

	def __str__(self):
		if not self.three: self._assign_three()
		self.rs.set_putative(self.three[0], self.three[1], self._score_five())	
		return str(self.rs)

files = ["../reads/timecourse/10_repA.bam", "../reads/timecourse/20_repA.bam", "../reads/timecourse/5_repA.bam",
		 "../reads/timecourse/total_repA.bam", "../reads/timecourse/10_repB.bam", "../reads/timecourse/20_repB.bam",
		 "../reads/timecourse/5_repB.bam", "../reads/timecourse/total_repB.bam", "../reads/timecourse/10_repC.bam",
		 "../reads/timecourse/20_repC.bam", "../reads/timecourse/5_repC.bam"]

genome_seq = load_genome(open(sys.argv[2], 'r'))
fp_pwm, tp_pwm= make_pwm(sys.argv[3], genome_seq)
f_min, f_max = get_min_score(fp_pwm), get_max_score(fp_pwm)
jxns = []
for line in open(sys.argv[1], 'r'):
	if line[0] == 'C': continue
	jxns += [RS(Entry(line), genome_seq)]
samfiles = [pysam.AlignmentFile(f) for f in files]

MIN_INTRON = 100
INDEL_LEN = 5

threes = set()
for line in open(sys.argv[3]):
	chrom, start, end, strand = line.strip().split()
	chrom = chrom[3:]
	three = int(start)+1 if strand == '-' else int(end)
	threes.add((chrom, strand, three))

# TODO:
# Make sure both splice sites have GT AG motifs
# check that sjr are going to the same place?

# Make sure sjr go to annotated 3'ss.


for samfile in samfiles:
	for jxn in jxns:
		for read in samfile.fetch(*jxn.get_coords()):
			if read.is_unmapped: continue
			if (jxn.rs.strand == '+') != (read.is_read1 == read.is_reverse): continue
			blocks = merge_blocks(read.get_blocks())
			for i, j in zip(blocks[:-1], blocks[1:]):
				begin, stop = i[1], j[0]
				if stop - begin < MIN_INTRON: continue
				b_m, s_m = genome_seq[jxn.rs.chrom][begin:begin+2], genome_seq[jxn.rs.chrom][stop-2:stop]
				if jxn.rs.strand == '+':
					if b_m != 'GT' or s_m != 'AG': continue
					five, three = begin, stop
				else:
					if revcomp(s_m) != 'GT' or revcomp(b_m) != 'AG': continue
					five, three = stop, begin

 			if (jxn.rs.chrom, jxn.rs.strand, three) in threes:
				jxn.add_sjr(five, three)

for jxn in jxns:
	print jxn

