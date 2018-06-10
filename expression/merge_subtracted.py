import sys

cur_chrom, cur_start, cur_end, cur_name, cur_count, cur_strand = '', '', '', '', 0, ''
for line in sys.stdin:
	chrom, start, end, name, count, strand = line.strip().split('\t')

	if cur_chrom and chrom == cur_chrom and 5 > abs(int(start) - int(cur_start)) and 5 > abs(int(end) - int(cur_end)) and strand == cur_strand:
		cur_count += int(count)

	else:
		if cur_chrom:
			print '\t'.join([cur_chrom, cur_start, cur_end, cur_name, str(cur_count), cur_strand])

		cur_chrom, cur_start, cur_end, cur_name, cur_count, cur_strand = (
			chrom, start, end, name, int(count), strand)

print '\t'.join([cur_chrom, cur_start, cur_end, cur_name, str(cur_count), cur_strand])
