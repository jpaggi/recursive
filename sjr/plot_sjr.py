import sys
import matplotlib.pyplot as plt
from core_pipeline.get_motifs import *
from core_pipeline.load_genome import *


def increased(peak, expression):
    WIDTH = 100
    before = sum(expression[peak - WIDTH: peak])
    after = sum(expression[peak:peak + WIDTH])

    return after > 2 * before
graveley = open('../data/graveley.bed', 'r')

grav = []
for line in graveley:
    chrom, start, end, name, score, strand = line.strip().split('\t')
    if strand == '+':
        grav += [(chrom, strand, int(end))]
    else:
        grav += [(chrom, strand, int(start))]

novel = open('../data/all_joined.bed', 'r')

nov = []
for line in novel:
    chrom, start, end, name, score, strand, seq1, seq2 = line.strip().split('\t')

    if seq1[-2:] != 'AG' or seq2[:2] != 'GT': continue
    if strand == '+':
        nov += [(chrom, strand, int(end))]
    else:
        nov += [(chrom, strand, int(start))]

genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

straddle = open('../data/pretty_big_straddle.bed', 'r')

strad = {}
for line in straddle:
    chrom, start, end, name, five, strand = line.strip().split('\t')

    if (chrom, strand, int(five) -1) in strad:
        strad[(chrom, strand, int(five) - 1)] += [int(end) if strand == '+' else int(start)]
    else:
        strad[(chrom, strand, int(five) - 1)] = [int(end) if strand == '+' else int(start)]
print strad

def color(p):
    if p > .9:
        return 'r'
    if p > .875: 
        return 'y'
    if p > .85:
        return 'k'
    if p > .8:
        return 'b'


data = open(sys.argv[1], 'r')
for line in data:

    chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]
    expression = [int(i) for i in line.strip().split('\t')[6].split(",")]
    start, end = int(start), int(end)
    
    seq = genome_seq[chrom][int(start):int(end)]
    if strand == "-":
        expression.reverse()
        seq = revcomp(seq)

    if True:
        plt.plot(expression)
        print chrom, start, end, strand, rs



        for gchrom, gstrand, grs in grav:
            if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
                if strand == '+':
                    pos = grs - int(start)
                elif strand == '-':
                    pos = int(end) - grs
                plt.axvline(pos, linewidth=4, color='g')

        for gchrom, gstrand, grs in nov:
            if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
                if strand == '+':
                    pos = grs - int(start)
                elif strand == '-':
                    pos = int(end) - grs
                plt.axvline(pos, linewidth=2, color='r')

        five = int(start) if strand == '+' else int(end)
        if (chrom, strand, five) in strad:
            for e in strad[(chrom, strand, five)]:
                if strand == '+':
                    pos = e - int(start)
                elif strand == '-':
                    pos = int(end) - e
                plt.axvline(pos + 1000, linewidth=2, color='k')

        for i in range(30, len(seq) - 30):
            motif = seq[i - len(tp_pwm): i + len(fp_pwm)]

            score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)

            if score > .8:
                plt.scatter([i], [1], linewidths = [(score - .8) * 100], c = color(score))

        plt.autoscale(tight = True)
        plt.show(block = False)
        a = raw_input("enter to continue")
        plt.close()
