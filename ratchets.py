import pysam
from get_anno_splices import get_jxns

samfile = pysam.AlignmentFile('../data/star_5_sjr.bam.bam', 'rb')
sjr = pysam.AlignmentFile('../data/ratchets.bam', 'wb', template = samfile)

ss = get_jxns('../data/anno.ss')
c = 0
for read in samfile.fetch():
	c += 1
	if not c % 10000: print c
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
			if (chrom, strand_str, five) in ss:
				if three not in ss[(chrom, strand_str, five)]:
					if strand and three < max(ss[(chrom, strand_str, five)]):
						sjr.write(read)
					elif not strand and three > min(ss[(chrom, strand_str, five)]):
						sjr.write(read)