import sys
import matplotlib.pyplot as plt
from get_motifs import *
from load_genome import *
from peaks import get_peaks

WINDOWS = 100

genome_seq = load_genome(open(sys.argv[1]))
fp_pwm, tp_pwm = make_pwm(sys.argv[2], genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)


def color(p):
    if p > .9:
        return 'r'
    if p > .875: 
        return 'y'
    if p > .85:
        return 'k'
    if p > .8:
        return 'b'

c = 0
y, n = 0,0
with open(sys.argv[3]) as fp:
    for line in fp:
        c += 1
        chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]
        expression = [int(float(i)) for i in line.strip().split('\t')[6].split(",")]
        mcmc = [float(i) for i in line.strip().split('\t')[7].split(",")]

        seq = genome_seq[chrom][int(start):int(end)]
        if strand == "-":
            expression.reverse()
            seq = revcomp(seq)

        print chrom, start, end, strand, rs

        peaks = get_peaks(mcmc, max(expression), 0.08)
        for peak in peaks:
            peak.plot(len(expression))

        x = range(0, len(expression), WINDOWS)
        if len(mcmc) != len(x): x = range(0, len(expression)+1, WINDOWS)
        plt.plot(x, mcmc, linewidth=6, color= 'k')
        plt.plot(expression)
        for i in range(30, len(seq) - 30):
            motif = seq[i - len(tp_pwm): i + len(fp_pwm)]
            score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)
            if score > .8:
                plt.scatter([i], [max(expression)], marker = 'o', c = color(score))

        plt.autoscale(tight = True)
        plt.show(block = False)
        a = raw_input("enter to continue")
        plt.close()
print y, n
