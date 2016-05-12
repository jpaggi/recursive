import sys
from intron_rs import IntronRS

for line in sys.stdin:
	intron = IntronRS(line)

	print '\t'.join(map(str, [intron.chrom, intron.rs, intron.rs + 1, '.', '.', intron.strand]))