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

novel = open('../data/manual_sjr_rs.bed', 'r')

sjr = []
for line in novel:
    chrom, start, end, name, score, strand = line.strip().split('\t') # , count, intron_count, seq1, seq2, score1, score2

    #if seq1[-2:] != 'AG' or seq2[:2] != 'GT': continue

    sjr += [(chrom, strand, int(start))]

genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

sawtooth = open('../data/high_confidence_sawtooth.bed', 'r')

saw = []
for line in sawtooth:
    chrom, start, end, name, score, strand= line.strip().split('\t')

    #if seq1[-2:] != 'AG' or seq2[:2] != 'GT': continue
    saw += [(chrom, strand, int(start))]

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
out = open(sys.argv[2], 'w')
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

        for gchrom, gstrand, grs in sjr:
            if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
                print 'here'
                if strand == '+':
                    pos = grs - int(start)
                elif strand == '-':
                    pos = int(end) - grs
                plt.axvline(pos, ymin = 0, ymax = .5, linewidth=2, color='r')

        for gchrom, gstrand, grs in saw:
            if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
                if strand == '+':
                    pos = grs - int(start)
                elif strand == '-':
                    pos = int(end) - grs
                plt.axvline(pos, ymin = .5, ymax = 1, linewidth=2, color='m')

        for i in range(30, len(seq) - 30):
            motif = seq[i - len(tp_pwm): i + len(fp_pwm)]

            score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)

            if score > .8:
 #               marker = '*' if increased(i, expression) else 'o'
                marker = 'o'
                plt.scatter([i], [1], marker = marker, linewidths = [(score - .8) * 100], c = color(score))

        plt.autoscale(tight = True)
        plt.show(block = False)
        a = raw_input("enter to continue")
        if a: out.write(line)
        plt.close()
