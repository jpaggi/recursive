import sys

cur_chrom, cur_start, cur_end, cur_strand, cur_align, cur_samples = '', '', '', '', set(), set()
for line in sys.stdin:
	chrom, start, end, sample, score, strand, align = line.strip().split('\t')

	if cur_chrom == chrom and cur_start == start and cur_end == end and cur_strand == strand:
		cur_align.add(align)
		cur_samples.add(sample)
	else:
		if cur_chrom:
			print '\t'.join([cur_chrom, cur_start, cur_end, ','.join(list(cur_samples)), str(len(cur_align)), cur_strand])

		cur_chrom, cur_start, cur_end, cur_strand, cur_align, cur_samples = chrom, start, end, strand, set([align]), set([sample])

print '\t'.join([cur_chrom, cur_start, cur_end, ','.join(list(cur_samples)), str(len(cur_align)), cur_strand])
