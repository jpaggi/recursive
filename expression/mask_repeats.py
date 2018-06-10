import sys
from repeats import Repeats

data = open(sys.argv[1], 'r')
repeats = Repeats(sys.argv[2])

for line in data:
    chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]
    expression = [int(i) for i in line.strip().split('\t')[6].split(",")]
    start, end = int(start), int(end)
    if strand == "-":
        expression.reverse()
    expression = repeats.mask(chrom, start, end, strand, expression)
    if strand == "-":
        expression.reverse()
    print '\t'.join([chrom, str(start), str(end), offsets, rs, strand, ','.join(map(str, expression))])
