import sys
from get_motifs import *
from load_genome import *
#import matplotlib.pyplot as plt
WINDOWS = 100
import numpy as np
from random import randrange
from peaks import get_peaks
from math import exp


genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

half = 500
k = 6 / float(half)
data = open(sys.argv[1], 'r')

graveley = open('../data/graveley.bed', 'r')

grav = {}
for line in graveley:
    chrom, start, end, name, score, strand = line.strip().split('\t')
    if strand == '+':
        grav[(chrom, strand, int(end))] = 0
    else:
        grav[(chrom, strand, int(start))] = 0

def logistic(x, k, x_0):
    if abs(x - x_0) > 2000: return 0
    return 1 / float(1 + exp(-k * (x - x_0)))

for line in data:
    chrom, start, end, offsets, rs, strand, expression, mcmc = line.strip().split('\t')
    expression = map(int, expression.split(","))
    mcmc = map(float, mcmc.split(','))
    start, end = int(start), int(end)
    
    seq = genome_seq[chrom][int(start):int(end)]
    if strand == "-":
        expression.reverse()
        seq = revcomp(seq)

    if sum(expression) / float(int(end) - int(start)) == 0: continue

    # for gchrom, gstrand, grs in grav:
    #     if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
    #         if strand == '+':
    #             pos = grs - int(start)
    #         elif strand == '-':
    #             pos = int(end) - grs
    #         # plt.axvline(pos, linewidth=4, color='g')

    scores = [0] * 30
    for pos in range(30, len(seq) - 30):
        motif = seq[pos - len(tp_pwm): pos + len(fp_pwm)]
        if motif[len(tp_pwm) - 2: len(tp_pwm)+2] != 'AGGT':
            score = 0
        else:
            score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)
        scores += [score]
    scores += [0] * 30

    peaks = get_peaks(mcmc, max(expression), .08)

    for peak in peaks:
        best = (-1, 0)

        # left
        x_0 = peak.start - 6 / k

        for i in range(peak.start):
            score = max(0, scores[i] - .8) * logistic(i, k, x_0)

            if score > best[1]:
                best = (i, score)

        # right
        x_0 = peak.end + 6 / k
        for i in range(peak.end, len(expression)):
            score = max(0, scores[i] - .8) * logistic(i, -k, x_0)

            if score > best[1]:
                best = (i, score)

        # center
        for i in range(peak.start, peak.end):
            score = max(0, scores[i] - .8) * logistic(i, -k, x_0)
            if score > best[1]:
                best = (i, score)


        if strand == '+':
            rs = start + best[0]
            if (chrom, strand, rs) in grav:
                grav[(chrom, strand, rs)] = best[1]
        else:
            rs = end - best[0]

            if (chrom, strand, rs) in grav:
                grav[(chrom, strand, rs)] = best[1]



    #     print best[1]
    #     plt.axvline(best[0], c = 'r', linewidth=2)

    # for peak in peaks:
    #     plt.axvline(peak.start, c = 'm')
    #     plt.axvline(peak.end, c = 'm')


    # for i in range(30, len(seq) - 30):
    #     motif = seq[i - len(tp_pwm): i + len(fp_pwm)]

    #     score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)

    #     if score > .8:
    #         plt.scatter([i], [max(expression)], linewidths = [(score - .8) * 100])
    # plt.plot(expression)
    # plt.autoscale(tight = True)
    # plt.show(block= False)
    # a = raw_input('enter to continue')
    # plt.close()

print grav.values(), len(filter(lambda x: x != 0, grav.values()))
