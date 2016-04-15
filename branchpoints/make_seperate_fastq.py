import sys
from fastq_iter import FastQIter
from separate import separate

reads1 = sys.argv[1]
reads2 = sys.argv[2]

out1 = open(sys.argv[3], 'w')
out2 = open(sys.argv[4], 'w')

for ID, seq1, seq2, qual1, qual2 in FastQIter(reads1, reads2):
	entry1, entry2 = separate(ID, seq1, seq2, qual1, qual2)

	out1.write(entry1)
	out2.write(entry2)
