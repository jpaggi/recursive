import pysam
import sys

samfile = pysam.AlignmentFile(sys.argv[1], "rb")
introns =  open(sys.argv[2], 'r')

for intron in introns:
    chrom, intron_start, intron_end, name1, name2, strand = intron.strip().split('\t')
    intron_start = int(intron_start)
    intron_end = int(intron_end)
    expression = [0] * (intron_end - intron_start)
    for read in samfile.fetch(chrom, intron_start, intron_end):
        for block in read.get_blocks():
            start = max(intron_start, block[0])
            end = min(intron_end, block[1])
            for i in xrange(start, end):
                expression[i - intron_start] += 1
    if sum(expression):
        print '\t'.join(map(str, [chrom, intron_start, intron_end, name1, name2, strand]) + [','.join(map(str, expression))])