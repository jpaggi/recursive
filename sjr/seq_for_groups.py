import sys
from subprocess import Popen, PIPE

SEQ_LENGTH = 30

introns = open(sys.argv[1], 'r')
genome_fasta = sys.argv[2]
out = open(sys.argv[3], 'w')

bedtools = Popen("bedtools getfasta -s -fi {} -bed - -fo stdout".format(genome_fasta), shell = True, stdin = PIPE, stdout = PIPE)

for intron in introns:
	chrom, start, end, sample, offsets, strand = intron.strip().split('\t')
	three = int(end) if strand == '+' else int(start)
	bedtools.stdin.write('\t'.join(map(str, [chrom, three - SEQ_LENGTH, three + SEQ_LENGTH, three, 0, strand])) + '\n')

introns.close()
bedtools.stdin.close()

bed = open(sys.argv[1], 'r')
fasta = bedtools.stdout

for line in bed:
	chrom, start, end, sample, offsets, strand = line.strip().split('\t')
	fasta_chrom, rest = fasta.readline().strip().split(':')
	fasta_chrom = fasta_chrom[1:]
	fasta_strand = rest[-2]
	rest = rest[:-3]
	fasta_start, fasta_end = rest.split('-')

	assert fasta_chrom == chrom
	if strand == '+':
		assert fasta_start == str(int(end) - SEQ_LENGTH)
	else:
		assert fasta_start == str(int(start) - SEQ_LENGTH)

	seq = fasta.readline().strip()

	out.write('\t'.join([chrom, start, end, sample, offsets, strand, seq[:SEQ_LENGTH], seq[SEQ_LENGTH:]]) + '\n')
