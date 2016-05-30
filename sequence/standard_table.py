"""
coords   site   graveley   (# of reads in junctions)   sawtooth   manual    motif_score

Other ideas...
Number of outgoing reads, exon body read counts, score based on increasing read count,
"null exon"
"""
from intron_exon import IntronRS
from load_genome import *
from get_motifs import *

sjr      = open('../data/rs_jxns/all_exons_counts.bed', 'r')
sawtooth = open('../data/mcmc_peak_calls/all_merged_masked.bed', 'r')
graveley = open('../data/graveley.bed', 'r')
exons    = open('../data/exons.bed')

seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta'))
fp_pwm, tp_pwm= make_pwm('../data/anno.ss', seq)
fp_min_score, fp_max_score = get_min_score(fp_pwm), get_max_score(fp_pwm)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

CSV = lambda x: ','.join(map(str, list(x)))

class Entry:
	def __init__(self, chrom, rs, strand):
		self.five = set()
		self.chrom = chrom
		self.rs = rs
		self.strand = strand
		self.motif = self._get_motif()
		self.three = 0
		self.grav = 0
		self.junc = [0] * 8
		self.sawtooth = 0
		self.sawtooth_prob = 0
		self.intron_sjr = 0
		self.motif_score = (score_motif(pwm, self.motif[5:5+28]) - min_score) / (max_score - min_score)
		self.manual = 0
		self.down_counts = [0] * 8
		self.body_counts = [0] * 4
		self.five_score = 0
		self.putative_three = 0
		self.putative_three_score = 0
		self.annotated_three = 0

	def add_auxillary_info(self, body_counts, down_counts, three):
		self.body_counts = body_counts
		self.down_counts = down_counts
		self.putative_three = three
		motif = seq[self.chrom][three - 8: three + 8]
		if self.strand == '-': motif = revcomp(motif)
		self.putative_three_score = (score_motif(fp_pwm, self.motif[-8:]) - fp_min_score) / (fp_max_score - fp_min_score)
		motif = seq[self.chrom][self.rs - 8: self.rs + 8]
		if self.strand == '-': motif = revcomp(motif)
		self.five_score = (score_motif(fp_pwm, motif[-8:]) - fp_min_score) / (fp_max_score - fp_min_score)

	def set_annotated_three(self, three):
		self.annotated_three = three

	def _get_motif(self):
		motif = seq[self.chrom][self.rs - 25:self.rs + 25]
		if self.strand == '-': motif = revcomp(motif)
		return motif[:37]

	def set_intron_sjr(self, count):
		self.intron_sjr = count

	def start(self):
		five = CSV(self.five) if self.five else '0'
		return five if self.strand =='+' else self.three

	def end(self):
		five = CSV(self.five) if self.five else '0'
		return self.three if self.strand =='+' else five

	def add_sjr(self, counts):
		self.junc = [i+j for i, j in zip(self.junc, counts)]

	def set_five(self, five):
		self.five.add(five)

	def set_three(self, three):
		self.three = three

	def add_sawtooth(self, score, best_prob):
		self.sawtooth = score
		self.sawtooth_prob = best_prob

	def add_graveley(self):
		self.grav = 1

	def add_manual(self, manual):
		self.manual = manual

	def _coords_str(self):
		return "{}:{}-{}:{}".format(self.chrom, self.start(), self.end(), self.strand)

	def recursive_index(self):
		return (sum(self.junc)+1) / float(self.intron_sjr+1)

	def __str__(self):
		return '\t'.join(map(str, [
			self._coords_str(),
			self.rs,
			self.grav,
			CSV(self.junc),
			self.sawtooth,
			self.sawtooth_prob,
			self.manual,
			self.intron_sjr,
			self.recursive_index(),
			self.motif,
			self.motif_score,
			CSV(self.down_counts),
			CSV(self.body_counts),
			self.five_score,
			self.putative_three_score,
			self.putative_three,
			self.annotated_three]))

entries = {}

for line in sjr:
	rs = IntronRS(line)
	key = (rs.chrom, rs.rs, rs.strand)
	if key not in entries:
		entries[key] = Entry(*key)
	entries[key].add_sjr(rs.counts)
	entries[key].set_five(rs.five())
	entries[key].add_auxillary_info(rs.body_counts, rs.down_counts, rs.put_three)

for line in sawtooth:
	chrom, rs, rs2, prob, score, strand = line.strip().split()

	key = (chrom, int(rs), strand)

	if key not in entries:
		entries[key] = Entry(*key)
	entries[key].add_sawtooth(score, prob)

for line in graveley:
	chrom, start, end, a, count, strand = line.strip().split()
	rs = int(end) if strand == '+' else int(start)
	five = int(start) if strand == '+' else int(end)
	key = (chrom, rs, strand)

	if key not in entries:
		entries[key] = Entry(*key)
	entries[key].add_graveley()
	entries[key].set_five(five)

for line in exons:
	chrom, start, end, pos, count, strand = line.strip().split()
	start, end = int(start), int(end)
	for entry in entries.values():
		if chrom == entry.chrom and strand == entry.strand:
			if strand == '+' and entry.rs == start:
				entry.set_annotated_three(end)
			elif strand == '-' and entry.rs == end:
				entry.set_annotated_three(start)

"""
assigns intron by:

1) highest number of splice jxn reads
2) longest intron

Overlooks introns shorted than 1000 bp (what ever is in intron file)
"""

# PRINT HEADER
print '\t'.join(['COORDS', 'RS', 'GRAV', 'SJR', 'SAW_SCORE', 'MCMC_PROB', 'MANUAL', 'SPANNING',
	             'RECURSIVE_INDEX' 'MOTIF', 'MOTIF_SCORE', 'DOWN', 'BODY', 'FIVE_SCORE', 'PUT_FIVE_SCORE',
	             'PUT_THREE', 'ANNO_THREE'])

introns = open('../data/all_merged.bed', 'r')

anno = {}
for line in introns:
	chrom, start, end, name, juncs, strand = line.strip().split('\t')[:6]
	start, end, juncs = int(start), int(end), int(juncs)
	if strand == '+':
		five, three = start, end
	else:
		five, three = end+1, start
	key = (chrom, five, strand)

	if key in anno:
		anno[key] += [(three, juncs)]
	else:
		anno[key] = [(three, juncs)]

for entry in entries.values():
	chrom, five, three, rs, strand = entry.chrom, entry.five, entry.three, entry.rs, entry.strand

	if not five:
		#assign five and three
		best = (-1, -1, -1)
		for key in anno:
			a_chrom, a_five, a_strand = key
			if a_chrom != chrom or a_strand != strand: continue
			for a_three, a_juncs in anno[key]:
				if a_juncs < best[2]: continue
				if strand == '+' and a_five + 5 < rs < a_three - 5:
					best = (a_five, a_three, a_juncs)
				elif strand == '-' and a_three + 5 < rs < a_five - 5:
					best = (a_five, a_three, a_juncs)

		if best[2] != -1:
			entry.set_five(best[0])
			entry.set_three(best[1])
			entry.set_intron_sjr(best[2])
			print entry
		else:
			assert False
	else:
		# assign three
		best = (0, -1)
		for f in five:
			key = (chrom, f, strand)
			if key in anno:
				for three, juncs in anno[key]:
					if juncs < best[1]: continue
					if strand == '+':
						if rs >= three + 5: continue
						if three > best[0] or juncs > best[1]: best = (three, juncs)
					if strand == '-':
						if rs + 5 <= three: continue
						if three < best[0] or juncs > best[1]: best = (three, juncs)
		if best[1] != -1:
			entry.set_three(best[0])
			entry.set_intron_sjr(best[1])
			print entry
		else:	
			# introns less than 1000 not present in expression file
			assert abs(f - rs) < 1000

