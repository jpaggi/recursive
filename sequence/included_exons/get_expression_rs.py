import pysam
from intron_rs import IntronRS
import sys

TIMES = ['5', '10', '20', 'total']
INDEX = lambda x: TIMES.index(x)

class RS:
	def __init__(self, intron_rs, three):
		self.rs = intron_rs
		self.three = int(three)

		self.body_count = [0] * len(TIMES)
		self.out_count  = [0] * len(TIMES) * 2

	def add_sjr(self, sample):
		self.out_count[INDEX(sample) * 2] += 1

	def add_per(self, sample):
		self.out_count[INDEX(sample) * 2 + 1] += 1

	def add_body(self, sample):
		self.body_count[INDEX(sample)] += 1

	def chrom_start_end_strand(self):
		if self.rs.strand == '+':
			return self.rs.chrom, self.rs.rs, self.three, self.rs.strand
		else:
			return self.rs.chrom, self.three, self.rs.rs, self.rs.strand

	def __str__(self):
		return '\t'.join(map(str,[self.rs, self.three, ','.join(map(str, self.body_count)), ','.join(map(str, self.out_count))]))

OVERHANG = 8
MIN_INTRON = 35
INDEL_LEN = 5

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

files = ["reads/timecourse/10_repA.bam", "reads/timecourse/20_repA.bam", "reads/timecourse/5_repA.bam",
		 "reads/timecourse/total_repA.bam", "reads/timecourse/10_repB.bam", "reads/timecourse/20_repB.bam",
		 "reads/timecourse/5_repB.bam", "reads/timecourse/total_repB.bam", "reads/timecourse/10_repC.bam",
		 "reads/timecourse/20_repC.bam", "reads/timecourse/5_repC.bam"]

# Create entries
entries = []
for line in open(sys.argv[1], 'r'):
	intron_rs = IntronRS('\t'.join(line.split('\t')[:-1]))
	entries.append(RS(intron_rs, line.strip().split('\t')[-1]))

for time in TIMES:
	for file in files:
		if time not in file: continue
		samfile = pysam.AlignmentFile(file)
		for entry in entries:
			chrom, start, end, strand = entry.chrom_start_end_strand()
			down = entry.three
			for read in samfile.fetch(chrom, int(start) - 300, int(end) + 300):
				# first body_count + out_count
				if read.is_unmapped: continue
				if (strand == '+') != (read.is_read1 == read.is_reverse): continue
				blocks = merge_blocks(read.get_blocks())
				for block in blocks:
					if start <= block[0] + OVERHANG and block[1] <= end + OVERHANG:
						entry.add_body(time)
						break
				for i, j in zip(blocks[:-1], blocks[1:]):
					begin, stop = i[1], j[0]
					if stop - begin < MIN_INTRON: continue
					if strand == '+':
						five, three = begin, stop
					else:
						five, three = stop, begin
					if five == down:
						entry.add_sjr(time)

				# now check if straddle read
				if not read.is_read2: continue
				if read.reference_id != read.next_reference_id: continue

				if strand == '+':
					inner_left  = blocks[-1][1]
					inner_right = read.pnext
				else:
					inner_left  = read.pnext
					inner_right = read.pos

				if inner_right - inner_left < 1000: continue
				if not (inner_left - OVERHANG < three < inner_right + OVERHANG): continue

				if (strand == '+') and end - inner_left < 300:
					entry.add_per(time)
				elif (not (strand == '+')) and inner_right - start < 300:
					entry.add_per(time)

for entry in entries:
	print entry
