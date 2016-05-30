
import sys
import pysam
import matplotlib.pyplot as plt

"""
Plots the coverage in transformed gene models.

Purely for debugging purposes.
"""

inside = True
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
introns =  open(sys.argv[2], 'r')

def merge_blocks(blocks):
    """
    Removes short gaps in read that are likely indels
    """
    i = 1
    while i < len(blocks) - 1:
        if blocks[i][0] - blocks[i-1][1] < 5:
            blocks[i] = (blocks[i][0], blocks[i+1][1])
            blocks.remove(blocks[i+1])
        else:
            i += 1
    return blocks

class Entry:
    def __init__(self, file):
        gene = file.readline()
        file.readline()
        file.readline()
        file.readline()
        exon1 = file.readline()
        exon2 = file.readline()

        self.chrom, a, event, self.start, self.end = gene.split('\t')[:5]
        h, a, g, b, self.exon1 = exon1.split('\t')[:5]
        h, a, c, self.exon2, b = exon2.split('\t')[:5]

        self.start = int(self.start)
        self.end = int(self.end)
        self.exon1 = int(self.exon1)
        self.exon2 = int(self.exon2)

while True:
    entry = Entry(introns)
    expression = [0] * (entry.end - entry.start)
    for read in samfile.fetch(entry.chrom, entry.start, entry.end):
        blocks = merge_blocks(read.get_blocks())
        # add expression for current read end
        for block in blocks:
            start = max(entry.start, block[0])
            end = min(entry.end, block[1])
            for i in xrange(start, end):
                expression[i - entry.start] += 1

    print entry.chrom
    print len(expression)
    plt.plot(expression)
    plt.axvline(entry.exon1 - 1000, linewidth = 2, c = 'r')
    plt.axvline(entry.exon2 - 1000, linewidth = 2, c = 'r')
    plt.show(block = False)

    a = raw_input('enter')

    plt.close()
