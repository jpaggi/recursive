import sys
import matplotlib.pyplot as plt
from load_genome import revcomp

class IntronRS:
	def __init__(self, line, seq = None):
		attributes = line.strip().split('\t')

		self.chrom   = attributes[0]
		self.start   = int(attributes[1])
		self.end     = int(attributes[2])
		self.rs      = int(attributes[3])
		self.counts  = map(int, attributes[4].split(','))
		self.strand  = attributes[5]
		self.intron  = int(attributes[6])

		if seq == None:
			self.motif = None
		else:
			self.motif = seq[self.chrom][self.rs-30:self.rs+30]
			if self.strand == '-': self.motif = revcomp(self.motif)

	def same(self, other):
		return (self.chrom  == other.chrom
			and self.start  == other.start
			and self.end    == other.end
			and self.strand == other.strand
			and self.rs     == other.rs)

	def merge(self, other):
		assert self.same(other)
		self.counts = [i+j for i, j in zip(self.counts, other.counts)]

	def __str__(self):
		return '\t'.join(map(str, [
				self.chrom, 
				self.start,
				self.end,
				self.rs,
				','.join(map(str, self.counts)),
				self.strand,
				self.intron,
			]))

	def length(self):
		return max(self.end - self.start, 0)

	def five(self):
		return self.start if self.strand == '+' else self.end

	def three(self):
		return self.start if self.strand == '-' else self.end

	def aggt(self):
		return self.motif[28:32] == 'AGGT'

	def motif_str(self, tp_len = 25, fp_len = 10):
		return self.seq1[30-tp_len:] + self.seq2[:30+fp_len]

	def motif_score(self, pwm, min_score, max_score):
		return (score_motif(self.motif_str(20, 8)) - min_score) / (max_score - min_score)

	def read_count(self):
		return sum(self.counts)

	def num_samples(self):
		return sum(map(lambda x: x > 0, self.counts))

	def expression(self):
		return self.intron

	def recursive_index(self):
		return self.read_count() / float(self.intron)

	def key(self):
		return (self.chrom, self.five(), self.strand)


