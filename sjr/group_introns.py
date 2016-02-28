import sys

if sys.argv[1] == '-':
	reads = sys.stdin
else:
	reads = open(sys.argv[1], 'r')
introns = open(sys.argv[2], 'w')

sample = sys.argv[3]


offsets = set()
cur_chrom, cur_start, cur_end, cur_strand = '', '', '', ''

for read in reads:
	chrom, start, end, name, blocks, strand = read.strip().split('\t')
	chrom = chrom[3:]
	blocks = eval(blocks)
	if cur_chrom and chrom == cur_chrom and start == cur_start and end == cur_end and strand == cur_strand:
		offsets.add((blocks[0][0], blocks[-1][1]))
	else:
		if cur_chrom:
			introns.write('\t'.join(map(str, [cur_chrom, cur_start, cur_end, sample, len(offsets), cur_strand])) + '\n')

		cur_chrom, cur_start, cur_end, cur_strand = chrom, start, end, strand
		offsets = set()
		offsets.add((blocks[0][0], blocks[-1][1]))

introns.write('\t'.join(map(str, [cur_chrom, cur_start, cur_end, sample, len(offsets), strand])) + '\n')
introns.close()