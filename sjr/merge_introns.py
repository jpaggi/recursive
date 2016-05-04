import sys

cur_chrom, cur_start, cur_end, cur_rs, cur_strand, cur_intron_count = '', '', '', [], '', ''
for line in sys.stdin:
	chrom, start, end, samples, rs, strand, count, intron_count, seq1, seq2, score1, score2 = line.strip().split('\t')

	if cur_chrom and chrom == cur_chrom and cur_start == start and cur_end == end and cur_strand == strand:
		cur_rs += [(int(rs), samples, count)]

	else:
		if cur_chrom and cur_rs:
			cur_rs = sorted(cur_rs)
			rs_str = ','.join(map(lambda x: str(x[0]), cur_rs))
			counts_sjr = ','.join(map(lambda x: str(x[2]), cur_rs))

			print '\t'.join([cur_chrom, cur_start, cur_end, counts_sjr, rs_str, cur_strand, cur_intron_count])


		cur_chrom, cur_start, cur_end, cur_rs, cur_strand, cur_intron_count = (
			chrom, start, end, [(int(rs), samples, count)], strand, intron_count)


cur_rs = sorted(cur_rs)
rs_str = ','.join(map(lambda x: str(x[0]), cur_rs))
counts_sjr = ','.join(map(lambda x: str(x[2]), cur_rs))
print '\t'.join([cur_chrom, cur_start, cur_end, counts_sjr, rs_str, cur_strand, cur_intron_count])
