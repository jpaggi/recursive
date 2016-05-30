import sys

reads1 = open(sys.argv[1], 'r')
reads2 = open(sys.argv[2], 'r')

discordant = open(sys.argv[3], 'w')

fastq1 = open(sys.argv[4], 'w')
fastq2 = open(sys.argv[5], 'w')

def process_aligns(aligns):
    out = []
    for align in aligns.split(' '):
        chrom, pos, cigar, strand = align.split(',')
        if chrom == '*':
            assert len(aligns.split(' ')) == 1
            return None
        out += [(chrom, pos, cigar, strand)]

    return out

for line1, line2 in zip(reads1, reads2):
    ID1, seq1, qual1, aligns1 = line1.strip().split('\t')
    ID2, seq2, qual2, aligns2 = line2.strip().split('\t')

    assert ID1 == ID2
    
    aligns1, aligns2 = process_aligns(aligns1), process_aligns(aligns2)

    # Neither end aligns ---> trash
    if aligns1 == aligns2 == None:
        continue

    # One end aligns ---> look for jxn read
    elif aligns1 == None:
        # Break read1
        #ID
        fastq1.write("@{}-{}\n".format(ID1, seq1))
        fastq2.write("@{}-{}\n".format(ID1, seq1))
        #SEQ
        fastq1.write(seq2 + '\n')
        fastq2.write(seq1[:15] + '\n')
        #BLANK
        fastq1.write('+\n')
        fastq2.write('+\n')
        #QUAL
        fastq1.write(qual2 + '\n')
        fastq2.write(qual1[:15] + '\n')


    elif aligns2 == None:
        # Break read2
        #ID
        fastq1.write("@{}-{}\n".format(ID1, seq2))
        fastq2.write("@{}-{}\n".format(ID1, seq2))
        #SEQ
        fastq1.write(seq2[:15] + '\n')
        fastq2.write(seq1 + '\n')
        #BLANK
        fastq1.write('+\n')
        fastq2.write('+\n')
        #QUAL
        fastq1.write(qual2[:15] + '\n')
        fastq2.write(qual1 + '\n')

    # Both ends align ---> check if properly discordant
    else:
        # Must be unique!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if len(aligns1) > 1 or len(aligns2) > 1: continue
        for chrom1, pos1, cigar1, strand1 in aligns1:
            for chrom2, pos2, cigar2, strand2 in aligns2:

                if chrom1 != chrom2: continue
                if strand1 == strand2: continue
                
                if strand1 == '-' and 0 < - int(pos1) + int(pos2) < 2000000:
                    discordant.write('\t'.join([chrom1, pos2, pos1, ID1, 'pair', strand1]) + '\n')
                
                if strand1 == '+' and 0 < - int(pos2) + int(pos1) < 2000000:
                    discordant.write('\t'.join([chrom1, pos1, pos2, ID1, 'pair', strand1]) + '\n')
