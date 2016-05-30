import sys
from sam_iter import SamIter
from load_genome import load_genome

def complement(char):
	if char == 'A': return 'T'
	if char == 'C': return 'G'
	if char == 'G': return 'C'
	if char == 'T': return 'A'
	return char

def revcomp(seq):
	return ''.join(map(complement, seq[::-1]))

def unique_pos(reads):
	out = []
	for test in reads:
		present = False
		for included in out:
			if test.rname == included.rname or test.pos == included.pos:
				present = True
		if not present:
			out.append(test)
	return out


def extend_from_end(genome_seq, broken_seq):
	N_ALLOWED_MISMATCHES = 2
	i = len(broken_seq) - 1
	mismatches = 0
	while mismatches < N_ALLOWED_MISMATCHES and i >= 0:
		if genome_seq[i] != broken_seq[i]:
			mismatches += 1
		i -= 1
	return i

def extend_from_start(genome_seq, broken_seq):
	N_ALLOWED_MISMATCHES = 2
	mismatches = 0
	j = 0
	while mismatches < N_ALLOWED_MISMATCHES and j < len(broken_seq):
		if genome_seq[j] != broken_seq[j]:
			mismatches += 1
		j += 1
	return j

def extend(read1, read2, genome):
	# Indexing not currently right!!!!!!!!!!

	if read1.is_reverse == read2.is_reverse: return False
	if read1.qname != read2.qname: return False
	if read1.rname in ['U', 'Uextra', 'dmel_mitochondrion_genome']: return False

	genome_seq = genome[read1.rname]
	broken_seq = revcomp(read1.qname.strip().split('-')[-1])

	INSERT_LENGTH = 400
	READ_LENGTH = len(broken_seq)
	SHORT_READ_LENGTH = min(len(read1.qual), len(read2.qual))
	MAX_INTRON_LENGTH = 2000000

	# read2 broken
	if len(read1.seq) == READ_LENGTH:
		if (not read1.is_reverse) and 0 < read1.pos - read2.pos < MAX_INTRON_LENGTH:
			# intron on plus strand
			i = extend_from_end(genome_seq[read2.pos+SHORT_READ_LENGTH-READ_LENGTH:read2.pos+SHORT_READ_LENGTH], broken_seq)
			if i > READ_LENGTH - SHORT_READ_LENGTH: return 'too short'
			if i < 1: return False#assert False, 'full length match'
			for shift in range(read1.pos + READ_LENGTH, read1.pos + READ_LENGTH + INSERT_LENGTH):
				j = extend_from_start(genome_seq[shift: shift+READ_LENGTH], broken_seq)
				if j >= i:
					for x in range(i, j+1):
						if broken_seq[x:x+2] == 'GT':
							print '1+', broken_seq[:x], broken_seq[x:], read1.rname, read1.pos, read2.pos+SHORT_READ_LENGTH-x-1, shift+x
		elif read1.is_reverse and 0 < read2.pos - read1.pos < MAX_INTRON_LENGTH:
			# intron on minus strand
			i = extend_from_end(revcomp(genome_seq[read2.pos:read2.pos+READ_LENGTH]), broken_seq)
			if i > READ_LENGTH - SHORT_READ_LENGTH: return 'too short'
			if i < 1: return False#assert False, 'full length match'
			for shift in range(read1.pos-INSERT_LENGTH, read1.pos):
				j = extend_from_start(revcomp(genome_seq[shift-READ_LENGTH: shift]), broken_seq)
				if j >= i:
					for x in range(i, j+1):
						if broken_seq[x:x+2] == 'GT':
							print '1-', broken_seq[:x], broken_seq[x:], read1.rname, read1.pos, read2.pos+READ_LENGTH-x, shift-x	

	# read1 broken
	elif len(read2.seq) == READ_LENGTH:
		if (not read2.is_reverse) and 0 < read2.pos - read1.pos < MAX_INTRON_LENGTH:
			# intron is on minus strand
			i = extend_from_end(genome_seq[read1.pos+SHORT_READ_LENGTH-READ_LENGTH:read1.pos+SHORT_READ_LENGTH], broken_seq)
			if i > READ_LENGTH - SHORT_READ_LENGTH: return 'too short'
			if i < 1: return False#assert False, 'full length match'
			for shift in range(read2.pos+READ_LENGTH, read2.pos+READ_LENGTH+INSERT_LENGTH):
				j = extend_from_start(genome_seq[shift: shift+READ_LENGTH], broken_seq)
				if j >= i:
					for x in range(i, j+1):
						if broken_seq[x-2:x] == 'AC':
							print '2-', broken_seq[:x], broken_seq[x:], read1.rname, read2.pos, read1.pos+SHORT_READ_LENGTH-x-3, shift+x
		elif read2.is_reverse and 0 < read1.pos - read2.pos < MAX_INTRON_LENGTH:
			# intron is on plus strand
			i = extend_from_end(revcomp(genome_seq[read1.pos:read1.pos+READ_LENGTH]), broken_seq)
			if i > READ_LENGTH - SHORT_READ_LENGTH: return 'too short'
			if i < 1: return False#assert False, 'full length match'
			for shift in range(read2.pos-INSERT_LENGTH, read2.pos):
				j = extend_from_start(revcomp(genome_seq[shift-READ_LENGTH: shift]), broken_seq)
				if j >= i:
					for x in range(i, j+1):
						if broken_seq[x-2:x] == 'AC':
							print '2+', broken_seq[:x], broken_seq[x:], read1.rname, read2.pos, read1.pos+READ_LENGTH-x, shift-x



def pair_reads(reads, genome):
	read1s = unique_pos(filter(lambda x: x.is_read1, reads))
	read2s = unique_pos(filter(lambda x: not x.is_read1, reads))

	aligns = []
	for read1 in read1s:
		for read2 in read2s:
			aligns.append(extend(read1, read2, genome))
	return aligns

			

genome = load_genome(open(sys.argv[2], 'r'))

cur_ID = ''
cur_reads = []
for read in SamIter(sys.argv[1]):
	if read.qname == cur_ID:
		cur_reads.append(read)
	else:
		pair_reads(cur_reads, genome)

		cur_ID = read.qname
		cur_reads = [read]

pair_reads(cur_reads, genome)
