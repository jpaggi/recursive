import sys
from get_motifs import *
from load_genome import *
WINDOWS = 100
from peaks import get_peaks
from math import exp


genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

half = 500
k = 6 / float(half)
data = open(sys.argv[1], 'r')

def logistic(x, k, x_0):
    if abs(x - x_0) > 2000: return 0
    return 1 / float(1 + exp(-k * (x - x_0)))

def call_peak(peak, k, scores):
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

    # arithmatic
    if best[1] > 0:
        if strand == '+':
            rs = start + best[0]
        else:
            rs = end - best[0]
        return (rs, best[1])
    return None

def call_peaks(peaks, k, scores):
    sites = []
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

        # arithmatic
        if best[1] > 0:
            if strand == '+':
                rs = start + best[0]
            else:
                rs = end - best[0]
            sites += [(rs, best[1])]
    return sites

for line in data:
    chrom, start, end, offsets, count, strand, expression, mcmc = line.strip().split('\t')
    expression = map(int, expression.split(","))
    mcmc = map(float, mcmc.split(','))
    start, end = int(start), int(end)
    
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


    # comment this line out to get real peaks!!!!!!!!
    # peaks = map(lambda x: x.random(len(expression)), peaks)

    for peak in peaks:
        call = call_peak(peak, k , scores)
        if call:
            rs, score = call
            print '\t'.join(map(str, [chrom, rs, rs+1, peak.best_prob, score, strand]))
