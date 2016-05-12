import sys
from get_motifs import *
from load_genome import *
import matplotlib.pyplot as plt
WINDOWS = 100
import numpy as np
from random import randrange
from peaks import get_peaks

def get_prob(pos, mcmc, z):
    window = abs(pos) / WINDOWS
    score = 0
    for i in range(-5, 10):
        try:
            score = max(score, mcmc[window + i])
        except IndexError: pass
    return score / float(z)

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

novel = open('../data/all_sjr_seq.bed', 'r')

nov = {}
for line in novel:
    chrom, start, end, name, score, strand, seq1, seq2, score1, score2 = line.strip().split('\t')

    if seq1[-2:] != 'AG' or seq2[:2] != 'GT': continue
    if strand == '+':
        nov[(chrom, strand, int(end))] = (seq1[-25:] + seq2[:10], score1, score2)
    else:
        nov[(chrom, strand, int(start))] = (seq1[-25:] + seq2[:10], score1, score2)

genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

straddle = open('../data/all_straddle_seq.bed', 'r')

strad = []
for line in straddle:
    chrom, start, end, name, score, strand, seq1, seq2 = line.strip().split('\t')[:8]

    if seq1[-2:] != 'AG' or seq2[:2] != 'GT': continue
    if strand == '+':
        strad += [(chrom, strand, int(end))]
    else:
        strad += [(chrom, strand, int(start))]

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

grav_probs = []
novel_probs = []
novel_random_probs = []
grav_random_probs = []

motif_probs = []
motif_random_probs = []

t_g, f, n, g, t_n, o, t_o = 0, 0, 0, 0, 0, 0, 0

n_w_p = 0
m = 0

bins = 21
step = 1 / float(bins - 1)

peaks_grav, peaks_novel, peaks_motif, motif_no_sjr = [0] * bins, [0] * bins, [0] * bins, [0] * bins
r_peaks_grav, r_peaks_novel, r_peaks_motif, r_motif_no_sjr = [0] * bins, [0] * bins, [0] * bins, [0] * bins

total_peaks = [0] * bins


for line in data:
    chrom, start, end, offsets, rs, strand, expression, mcmc = line.strip().split('\t')
    expression = map(int, expression.split(","))
    mcmc = map(float, mcmc.split(','))
    start, end = int(start), int(end)

    #if sum(expression) / len(expression) < 5: continue
    
    seq = genome_seq[chrom][int(start):int(end)]
    if strand == "-":
        expression.reverse()
        seq = revcomp(seq)

    if sum(expression) / float(int(end) - int(start)) == 0: continue

    peaks = []
    for thresh in [x * step for x in range(bins)]:
        peaks += [get_peaks(mcmc, max(expression), thresh)]
    random_peaks = []
    for peak in peaks:
        random_peaks += [map(lambda x: x.random(len(expression)), peak)]

    for i in range(len(peaks)):
        total_peaks[i] += len(peaks[i])
    

    for gchrom, gstrand, grs in grav:
        if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
            if strand == '+':
                pos = grs - int(start)
            elif strand == '-':
                pos = int(end) - grs
            grav_probs += [get_prob(pos, mcmc, max(expression))]
            grav_random_probs += [get_prob(randrange(100, len(expression)-100), mcmc, max(expression))]
            g += 1
            for i in range(len(peaks)):
                for peak in peaks[i]:
                    if peak.close(pos, 1000):
                        peaks_grav[i] += 1
                for peak in random_peaks[i]:
                    if peak.close(pos, 1000):
                        r_peaks_grav[i] += 1


    for gchrom, gstrand, grs in nov:
        if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
            if (gchrom, gstrand, grs) in grav:
                continue
            n += 1
            if strand == '+':
                pos = grs - int(start)
            elif strand == '-':
                pos = int(end) - grs
            for i in range(len(peaks)):
                for peak in peaks[i]:
                    if peak.close(pos, 1000):
                        peaks_novel[i] += 1
                for peak in random_peaks[i]:
                    if peak.close(pos, 1000):
                        r_peaks_novel[i] += 1
            novel_probs += [get_prob(pos, mcmc, max(expression))]
            novel_random_probs += [get_prob(randrange(100, len(expression)-100), mcmc, max(expression))]
            if get_prob(pos, mcmc, max(expression)) > .1:
                print nov[(gchrom, gstrand, grs)]
                n_w_p += 1

    # for gchrom, gstrand, grs in strad:
    #     if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
    #         if strand == '+':
    #             pos = grs - int(start)
    #         elif strand == '-':
    #             pos = int(end) - grs

    for pos in range(30, len(seq) - 30):
        motif = seq[pos - len(tp_pwm): pos + len(fp_pwm)]

        score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)
        if score > .85:
            m += 1
            for i in range(len(peaks)):
                for peak in peaks[i]:
                    if peak.close(pos, 1000):
                        peaks_motif[i] += 1
                for peak in random_peaks[i]:
                    if peak.close(pos, 1000):
                        r_peaks_motif[i] += 1
            motif_probs += [get_prob(pos, mcmc, max(expression))]
            motif_random_probs += [get_prob(randrange(100, len(expression)-100), mcmc, max(expression))]

def get_cdf(probs):
    y_vals = []
    x_vals = [x * 0.01 for x in range(100)]
    for x in x_vals:
        y_vals += [len(filter(lambda y: y >= x, probs))]
    return x_vals, y_vals

def plot_cdf(probs, random_probs, title, path):
    plt.title(title)
    plt.xlabel('score cutoff')
    plt.ylabel('Number of Events Detected')

    x_vals, y_vals = get_cdf(probs)
    c_x_vals, c_y_vals = get_cdf(random_probs)
    
    plt.plot(c_x_vals, c_y_vals, color = 'k', linewidth = 2)
    plt.plot(x_vals, y_vals, color = 'g', linewidth = 2)
    plt.savefig(path)
    plt.close()

def plot_mcmc(peaks, r_peaks):
    x_vals = [x * step for x in range(bins)]
    plt.plot(x_vals, r_peaks, color = 'k', linewidth = 2)
    plt.plot(x_vals, peaks, color = 'g', linewidth = 2)
    plt.show()



# size, Size = 'short', 'Short'
# plot_cdf(grav_probs, grav_random_probs, "Graveley {} Introns".format(Size), "../data/mcmc_evaluate/{}_graveley.png".format(size))
# plot_cdf(novel_probs, novel_random_probs, "Novel {} Introns".format(Size), "../data/mcmc_evaluate/{}_novel.png".format(size))
# plot_cdf(motif_probs, motif_random_probs, "Motif {} Introns".format(Size), "../data/mcmc_evaluate/{}_motif.png".format(size))

r_all_peaks = [i+j for i, j in zip(r_peaks_grav, r_peaks_novel)]
all_peaks = [i+j for i, j in zip(peaks_grav, peaks_novel)]

x_vals = [x * step for x in range(bins)]
plt.title('Short introns')
plt.ylabel('Number of Peaks')
plt.xlabel('Peak cutoff')
plt.plot(x_vals, r_all_peaks, color = 'k', linewidth = 2, label = 'Random for SJR')
plt.plot(x_vals, r_peaks_motif, color = 'y', linewidth = 2, label = 'Random for Motifs')
plt.plot(x_vals, peaks_grav, color = 'g', linewidth = 2, label  = 'Graveley')
plt.plot(x_vals, all_peaks, color = 'r', linewidth = 2, label = 'ALL SJRs and Graveley')
plt.plot(x_vals, peaks_motif, color = 'm', linewidth = 2, label = 'Motifs scoring > .8')
plt.plot(x_vals, total_peaks, linewidth= 2, label = 'All peaks')
plt.ylim([0, 500])
plt.legend()
plt.savefig('../data/mcmc_evaluate/short_recall.png')
plt.show()

print 'w/ p', len(filter(lambda x: x > .1, grav_probs))
print 'grav hits:', peaks_grav, r_peaks_grav
print 'graveley total', g
print

plot_mcmc(peaks_grav, r_peaks_grav)

print 'novel hits:', peaks_novel, r_peaks_novel
print 'w/ p', len(filter(lambda x: x > .1, novel_probs))
print 'novel total', n
print

plot_mcmc(peaks_novel, r_peaks_novel)

print 'motif hits:', peaks_motif, r_peaks_motif
print 'w/ p', len(filter(lambda x: x > .1, motif_probs))
print 'motifs total', m
print

plot_mcmc(peaks_motif, r_peaks_motif)

print 'total peaks:', total_peaks

