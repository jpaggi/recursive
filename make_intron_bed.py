import pysam
from get_anno_splices import get_three_prime_jxns

samfile = pysam.AlignmentFile('../data/star_5_3_prime.bam', 'rb')

ss = get_three_prime_jxns('../data/anno.ss')
c = 0
for read in samfile.fetch():
	c += 1
	blocks = read.get_blocks()
	if len(blocks) > 1:
		i = 1
		while i < len(blocks) - 1:
			if blocks[i][0] - blocks[i-1][1] < 5:
				blocks[i] = (blocks[i][0], blocks[i+1][1])
				blocks.remove(blocks[i+1])
			else:
				i += 1
		strand = (read.is_read1 == read.is_reverse)
		for i in xrange(len(blocks) - 1):
			if strand:
				five, three = blocks[i][1] - 1, blocks[i+1][0]
			else:
				three, five = blocks[i][1] - 1, blocks[i+1][0]
			chrom = 'chr' + samfile.getrname(read.reference_id)
			strand_str = '+' if strand else '-'
			if (chrom, strand_str, three) in ss:
				if five not in ss[(chrom, strand_str, three)]:
					if strand and five > min(ss[(chrom, strand_str, three)]):
						print '\t'.join(map(str, [chrom[3:], five, three, strand_str]))
					elif not strand and five < max(ss[(chrom, strand_str, three)]):
						print '\t'.join(map(str,[chrom[3:], three, five, strand_str]))