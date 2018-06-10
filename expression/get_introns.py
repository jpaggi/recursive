"""
This module takes as input a gtf file.

And prints introns in bed format to stdout.
* Uses half-open 0-indexing *

(Introns are defined as being regions between neightboring exons.)
"""

import sys

class Exon:
    def __init__(self, line):
        sep = line.split('\t')
        self.chrom = sep[0]
        self.source = sep[1]
        self.feature = sep[2]
        self.start = int(sep[3])
        self.end = int(sep[4])
        self.score = sep[5]
        self.strand = sep[6]
        self.transcript_name, self.gene_name = self.get_gene_name(sep[8])

    def get_gene_name(self, attributes):
        attributes = {attribute.split('"')[0].strip(): attribute.split('"')[1].strip() for attribute in attributes.split(';')[:-1]}
        return attributes['transcript_id'], attributes['gene_id']

# asserts that exons don't overlap and prints regions between them as introns
def print_transcript(transcript_name, transcript_exons):
    transcript_exons = sorted(transcript_exons, key = lambda x: x.start)
    for i in range(1, len(transcript_exons)):
        assert(transcript_exons[i- 1].end < transcript_exons[i].start)
        print '\t'.join([transcript_exons[i].chrom, str(transcript_exons[i - 1].end), str(transcript_exons[i].start - 1), transcript_exons[i].transcript_name, 
                         transcript_exons[i].gene_name, transcript_exons[i].strand])

processed = set()
transcript_name = ''
transcript_exons = []
with open(sys.argv[1]) as fp:
    for line in fp:
        if line.split('\t')[2] != 'exon': continue
        cur_exon = Exon(line)
        assert not cur_exon.transcript_name in processed

        if cur_exon.transcript_name != transcript_name:
            processed.add(transcript_name)
            print_transcript(transcript_name, transcript_exons)
            transcript_exons = [cur_exon]
            transcript_name = cur_exon.transcript_name
        else:
            transcript_exons += [cur_exon]

print_transcript(transcript_name, transcript_exons)
