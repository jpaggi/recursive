"""
This module takes as input only the exon lines of a gtf file.

*** file must be sorted by transcript name ***

I.e. all exons in a given transcript must be next to each other.
[This is hard to check so it is not checked in this module...]

Prints introns in bed format to stdout.

Introns are defined as being regions between neightboring exons.
"""

import sys

class Exon:
    def __init__(self, line):
        sep = line.split('\t')
        self.chrom = sep[0][3:]
        self.source = sep[1]
        self.feature = sep[2]
        self.start = int(sep[3])
        self.end = int(sep[4])
        self.score = sep[5]
        self.strand = sep[6]
        self.transcript_name, self.gene_name = self.get_gene_name(sep[8])

    def get_gene_name(self, attributes):
        return attributes.split('"')[1], attributes.split('"')[3]

# asserts that exons don't overlap and prints regions between them as introns
def print_transcript(transcript_name, transcript_exons):
    transcript_exons = sorted(transcript_exons, key = lambda x: x.start)
    for i in range(1, len(transcript_exons)):
        assert(transcript_exons[i- 1].end < transcript_exons[i].start)
        if transcript_exons[i].start - 2 - transcript_exons[i - 1].end < 1000: continue
        print '\t'.join([transcript_exons[i].chrom, str(transcript_exons[i - 1].end), str(transcript_exons[i].start - 2), transcript_exons[i].transcript_name, 
                         transcript_exons[i].gene_name, transcript_exons[i].strand])

#if transcript not contiguous in input will miss introns!
exons = sys.stdin
transcript_name = ''
transcript_exons = []
for line in exons:
    if line.split('\t')[2] != 'exon': continue
    cur_exon = Exon(line)
    if cur_exon.transcript_name != transcript_name:
        print_transcript(transcript_name, transcript_exons)
        transcript_exons = [cur_exon]
        transcript_name = cur_exon.transcript_name
    else:
        transcript_exons += [cur_exon]

print_transcript(transcript_name, transcript_exons)
