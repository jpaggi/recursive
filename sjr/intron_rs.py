import sys
import matplotlib.pyplot as plt

class IntronRS:
	def __init__(self, line):
		attributes = line.strip().split('\t')

		self.chrom   = attributes[0]
		self.start   = int(attributes[1])
		self.end     = int(attributes[2])
		self.samples = set(attributes[3].split(','))
		self.rs      = int(attributes[4])
		self.strand  = attributes[5]     
		self.sjr     = int(attributes[6])
		self.intron  = int(attributes[7])
		self.seq1    = attributes[8]
		self.seq2    = attributes[9]
		self.score1  = float(attributes[10])
		self.score2  = float(attributes[11])

	def same(self, other):
		return (self.chrom  == other.chrom
			and self.start  == other.start
			and self.end    == other.end
			and self.strand == other.strand
			and self.rs     == other.rs)

	def merge(self, other):
		assert self.same(other)
		self.sjr += other.sjr
		self.samples |= other.samples

	def __str__(self):
		return '\t'.join(map(str, [
				self.chrom, 
				self.start,
				self.end,
				','.join(self.samples),
				self.rs,
				self.strand,
				self.sjr,
				self.intron,
				self.seq1,
				self.seq2,
				self.score1,
				self.score2
			]))

	def length(self):
		return max(self.end - self.start, 0)

	def good_motifs(self):
		return self.score1 > .8 and self.score2 > .8

	def great_motifs(self):
		return self.score1 > .9 and self.score2 > .9

	def motifs_above_thresh(self, thresh):
		return self.score1 > thresh and self.score2 > thresh

	def motif_strengths(self):
		return self.score1, self.score2

	def aggt(self):
		return self.seq1[-2:] == 'AG' and self.seq2[:2] == 'GT'

	def motif_str(self, tp_len = 25, fp_len = 10):
		return self.seq1[-tp_len:] + self.seq2[:fp_len]

	def read_count(self):
		return self.sjr

	def num_samples(self):
		return len(self.samples)

	def expression(self):
		return self.intron

	def recursive_index(self):
		return self.sjr / float(self.intron)
