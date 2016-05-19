import sys
from load_genome import *

genome_seq = load_genome(open(sys.argv[2], 'r'))

exons = open(sys.argv[1], 'r')

out_motif = open(sys.argv[3], 'w')
out_no_motif = open(sys.argv[4], 'w')

for line in exons:
	chrom, start, end, position, score, strand = line.strip().split()[:6]

	try:
		if strand == '+':
			motif = genome_seq[chrom][int(start)-2: int(start)+2]
		else:
			motif = revcomp(genome_seq[chrom][int(end)-2: int(end)+2])
	except KeyError:
		continue

	if motif == 'AGGT':
		out_motif.write(line)
	else:
		out_no_motif.write(line)
