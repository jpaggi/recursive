import sys

reads1 = open(sys.argv[1], 'r')
reads2 = open(sys.argv[2], 'r')

out1 = open(sys.argv[3] + '1.fq', 'w')
out2 = open(sys.argv[3] + '2.fq', 'w')

def complement(char):
	if 'A': return 'T'
	if 'C': return 'G'
	if 'G': return 'C'
	if 'T': return 'A'
	return char

def get_read(reads, rev = False):
	ID = reads.readline()
	seq = reads.readline()
	ID2 = reads.readline()
	quality = reads.readline()
	if rev:
		seq = map(complement, seq[::-1])
		quality = quality[::-1]
	return ID, seq, quality




def get_pair(reads1, reads2):
	read1 = reads1.readline()
	read2 = reads2.readline()
