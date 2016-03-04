import sys
CONTEXT = 46

fasta = open(sys.argv[1], 'r')

def get_next(fasta):
	name = fasta.readline().strip()[1:]
	if not name:
		return '', '', '', '', ''
	seq = fasta.readline().strip()
	chrom, rest = name.split(':')
	strand = rest[-2]
	rest = rest[:-3]
	start, end = map(int, rest.split('-'))

	return chrom, start, end, strand, seq


while True:
	chrom, start, end, strand, seq = get_next(fasta)
	if not chrom: break

	for i in xrange(CONTEXT, len(seq) - CONTEXT - 3):
		if seq[i: i + 4] in ['AGGT', 'AGGC']:
			if strand == '+':
				five = start + CONTEXT - 1
				three = start + i + 2
			else:
				five = end - CONTEXT
				three = end - i - 3

			print ">{}:{}-{}({})".format(chrom, five, three, strand)
			print seq[:CONTEXT] + seq[i + 2: i + 2 + CONTEXT]
