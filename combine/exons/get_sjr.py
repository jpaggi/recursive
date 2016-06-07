import pysam
import sys

exons = open(sys.argv[1], 'r')
samfile = pysam.AlignmentFile(sys.argv[2])

OVERHANG = 8
MIN_INTRON = 35
INDEL_LEN = 5

def merge_blocks(blocks):
    """
    Removes short gaps in read that are likely indels
    """
    i = 1
    while i < len(blocks) - 1:
        if blocks[i][0] - blocks[i-1][1] < INDEL_LEN:
            blocks[i] = (blocks[i][0], blocks[i+1][1])
            blocks.remove(blocks[i+1])
        else:
            i += 1
    return blocks

for line in exons:
	chrom, start, end, position, score, strand = line.strip().split()
	body_count, up_count, down_count = 0, 0, 0
	if strand == '+':
		up, down = int(start), int(end)
	else:
		down, up = int(start), int(end)
	for read in samfile.fetch(chrom, int(start), int(end)):
		if read.is_unmapped: continue
		if (strand == '+') != (read.is_read1 == read.is_reverse): continue
		blocks = merge_blocks(read.get_blocks())
		for block in blocks:
			if int(start) <= block[0] + OVERHANG and block[1] <= int(end) + OVERHANG:
				body_count += 1
				break
		for i, j in zip(blocks[:-1], blocks[1:]):
			begin, stop = i[1], j[0]
			if stop - begin < MIN_INTRON: continue
			if strand == '+':
				five, three = begin, stop
			else:
				five, three = stop, begin
			if five == down:
				down_count += 1
			if three == up:
				up_count += 1
	print '\t'.join(map(str, [chrom, start, end, position, body_count, strand, up_count, down_count]))
