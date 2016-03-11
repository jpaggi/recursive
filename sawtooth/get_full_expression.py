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
        if (read.is_read2 and not read.is_reverse):
            expression[min(read.get_blocks()[-1][-1] - intron_start, len(expression) - 1)] += 1
        elif read.is_read1 and read.is_reverse:
            expression[max(read.get_blocks()[0][0] - intron_start, 0)] += 1
                
    if sum(expression):
        print '\t'.join(map(str, [chrom, intron_start, intron_end, name1, name2, strand]) + [','.join(map(str, expression))])