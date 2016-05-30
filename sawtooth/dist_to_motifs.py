import sys
from get_motifs import *
from load_genome import *
import matplotlib.pyplot as plt
WINDOWS = 100
import numpy as np
from random import randrange
from peaks import get_peaks

genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)


data = open(sys.argv[1], 'r')


dists = range(0, 5000, 200)

real_dists = []

fake_dists = []

peak_scores = [[] for dist in dists]
r_peak_scores = [[] for dist in dists]

peak_footprint = 0
intron_footprint = 0
m = 0
total_peaks = 0

MOTIF_CUTOFF = .85

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
    r_peaks = map(lambda x: x.random(len(seq)), peaks)

    intron_footprint += len(expression)

    for peak in peaks:
        peak_footprint += peak.length()
        total_peaks += 1

    for peak in peaks:
        if max(scores[peak.start: peak.end]) > MOTIF_CUTOFF:
            real_dists += [0]
        else:
            left, right = float('inf'), float('inf')
            for i in range(0, peak.start)[::-1]:
                if scores[i] > MOTIF_CUTOFF: 
                    left = i - peak.start
                    break
            for i in range(peak.end, len(seq)):
                if scores[i] > MOTIF_CUTOFF:
                    right = i - peak.end
                    break
            if abs(left) < abs(right):
                dist = left
            else:
                dist = right
            if dist == float('inf'):
                real_dists += [50000]
            else:
                real_dists += [dist]
        for i, dist in enumerate(dists):
            peak_scores[i] += [max(scores[max(0, peak.start - dist): min(peak.end + dist, len(seq))])]


    for peak in r_peaks:
        if max(scores[peak.start: peak.end]) > MOTIF_CUTOFF:
            fake_dists += [0]
        else:
            left, right = float('inf'), float('inf')
            for i in range(0, peak.start)[::-1]:
                if scores[i] > MOTIF_CUTOFF: 
                    left = i - peak.start
                    break
            for i in range(peak.end, len(seq)):
                if scores[i] > MOTIF_CUTOFF:
                    right = i - peak.end
                    break
            if abs(left) < abs(right):
                dist = left
            else:
                dist = right
            if dist == float('inf'):
                fake_dists += [50000]
            else:
                fake_dists += [dist]
        for i, dist in enumerate(dists):
            r_peak_scores[i] += [max(scores[max(0, peak.start - dist): min(peak.end + dist, len(seq))])]

    for pos in range(30, len(seq) - 30):
        motif = seq[pos - len(tp_pwm): pos + len(fp_pwm)]

        score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)
        if score > MOTIF_CUTOFF: m += 1


r_file = open('../data/long_peak_dists/random_dists.csv', 'w')

r_file.write(','.join(map(str, fake_dists)))

real_file = open('../data/long_peak_dists/real_dists.csv', 'w')

real_file.write(','.join(map(str, real_dists)))



r_score_file = open('../data/long_peak_dists/random_motif_scores.csv', 'w')

r_score_file.write(','.join(map(str, r_peak_scores[0])))

real_score_file = open('../data/long_peak_dists/real_motif_scores.csv', 'w')

real_score_file.write(','.join(map(str, peak_scores[0])))

counts = []
for i in range(0, 10000, 100):
    counts += [len(filter(lambda x: abs(x) <= i, real_dists))]
r_counts = []
for i in range(0, 10000, 100):
    r_counts += [m * (peak_footprint + total_peaks * i * 2) / float(intron_footprint)]
f_counts = []
for i in range(0, 10000, 100):
    f_counts += [len(filter(lambda x: abs(x) <= i, fake_dists))]


print 'real', len(filter(lambda x: x == 0, real_dists))
print 'fake', len(filter(lambda x: x == 0, fake_dists))
plt.plot(range(0, 10000, 100), counts)
plt.plot(range(0, 10000, 100), r_counts)
plt.plot(range(0, 10000, 100), f_counts)

plt.savefig('../data/long_peak_dists/dist_to_motif_cdf.png')
plt.show()

scores = []
for i in range(len(peak_scores)):
    scores += [[peak_scores[i]]]
    scores += [[r_peak_scores[i]]]
#scores = [[i]+[j] for i,j in zip(peak_scores, r_peak_scores)]
plt.boxplot(scores, positions = range(len(scores)))#, labels = dists)
#plt.boxplot(peak_scores, positions = range(0, 2 * len(peak_scores), 2), labels = dists)
#plt.boxplot(r_peak_scores, positions = range(1, 2 * len(peak_scores) + 1, 2), labels = dists)
plt.savefig('../data/long_peak_dists/best_motif_in_range.png')
plt.show()
plt.hist(peak_scores[0])
plt.show()
plt.hist(r_peak_scores[0])
plt.show()

