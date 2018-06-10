from load_genome import revcomp

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