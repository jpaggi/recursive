import sys
import matplotlib.pyplot as plt
WINDOWS = 100

from get_motifs import *
from load_genome import *


genome_seq = load_genome(open('../data/downloaded/dmel-all-chromosome-r5.57.fasta', 'r'))
fp_pwm, tp_pwm = make_pwm('../data/anno.ss', genome_seq)
pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)


def increased(peak, expression):
    WIDTH = 100
    before = sum(expression[peak - WIDTH: peak])
    after = sum(expression[peak:peak + WIDTH])

    return after > 2 * before

class Peak:
    def __init__(self, start, end):
        self.start = (start-1) * WINDOWS
        self.end = end * WINDOWS

    def inside(self, pos):
        return self.start <= pos <= self.end

    def close(self, pos, dist):
        return self.start - dist <= pos <= self.end + dist

    def length(self):
        return self.end - self.start

    def random(self, length):
        start = randrange(100, length)

    def plot(self, length):
        start_p = self.start / float(length)
        end_p = self.end / float(length)
        plt.axvline(self.start, c = 'm')
        plt.axvline(self.end, c = 'm')

    def mergable(self, end):
        end = end * WINDOWS
        return end < self.end + 500

    def merge(self, end):
        self.end = end * WINDOWS

def get_peaks(probs, z):
    """
    Given probs in reduced space
    Return Peaks in genome space
    """
    peaks = []
    peak_start = None
    THRESH = .1
    for i in xrange(len(probs)):
        p = probs[i] / float(z)
        # make new peak
        if peak_start != None and p < THRESH:
            if peaks and peaks[-1].mergable(i):
                peaks[-1].merge(i)
            else:
                peaks += [Peak(peak_start, i)]
            peak_start = None
        elif peak_start == None and p >= THRESH:
            peak_start = i
    return peaks


def get_probs(samples, expression):
    length = (len(expression) / WINDOWS) + 1
    counts = [0] * length
    for sample in samples:
        for ratchet in sample:
            counts[ratchet] += 1
    max_expression = max(expression)
    for i in xrange(len(counts)):
        counts[i] = (counts[i] * max_expression) / float(len(samples))
    return counts

def get_prob(pos, start, mcmc, z):
    window = abs(pos - start) / WINDOWS
    score = 0
    for i in range(-5, 10):
        try:
            score = max(score, mcmc[window + i])
        except IndexError: pass
    return score / float(z)

def color(p):
    if p > .9:
        return 'r'
    if p > .875: 
        return 'y'
    if p > .85:
        return 'k'
    if p > .8:
        return 'b'

graveley = open('../data/graveley.bed', 'r')

grav = []
for line in graveley:
    chrom, start, end, name, score, strand = line.strip().split('\t')
    if strand == '+':
        grav += [(chrom, strand, int(end))]
    else:
        grav += [(chrom, strand, int(start))]

novel = open('../data/all_sjr_seq.bed', 'r')

nov = []
for line in novel:
    chrom, start, end, name, score, strand, seq1, seq2 = line.strip().split('\t')[:8]

    if seq1[-2:] != 'AG' or seq2[:2] != 'GT': continue
    if strand == '+':
        nov += [(chrom, strand, int(end))]
    else:
        nov += [(chrom, strand, int(start))]


straddle = open('../data/all_straddle_seq.bed', 'r')

strad = []
for line in straddle:
    chrom, start, end, name, five, strand, seq1, seq2 = line.strip().split('\t')[:8]

    if seq1[-2:] != 'AG' or seq2[:2] != 'GT': continue
    if strand == '+':
        strad += [(chrom, strand, int(end))]
    else:
        strad += [(chrom, strand, int(start))]

data = open(sys.argv[1], 'r')
c = 0
y, n = 0,0
for line in data:
    c += 1
    #print c
    #if c < 200: continue 
    chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]
    expression = [int(float(i)) for i in line.strip().split('\t')[6].split(",")]
    mcmc = [float(i) for i in line.strip().split('\t')[7].split(",")]


    seq = genome_seq[chrom][int(start):int(end)]
    if strand == "-":
        expression.reverse()
        seq = revcomp(seq)

    peaks = get_peaks(mcmc, max(expression))

    # for peak in peaks:
    #     if increased(peak * WINDOWS, expression):
    #         plt.axvline(peak * WINDOWS, linewidth=3, color='r')
    #     else:
    #         plt.axvline(peak * WINDOWS, linewidth=3, color='m')

#     missed = False
#     for gchrom, gstrand, grs in grav:
#         if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
#             if strand == '+':
#                 p = get_prob(grs, int(start), mcmc, max(expression))
#                 pos = grs - int(start)
# #                plt.axvline(int(grs) - int(start), linewidth=4, color='g')
#             elif strand == '-':
#                 p = get_prob(grs, int(end), mcmc, max(expression))
#                 pos = int(end) - grs
# #                plt.axvline(int(end) - int(grs), linewidth=4, color='g')
#             if p > .1:
#                 y += 1
#             else:
#                 n += 1
#                 missed = True
            # inside = False
            # for peak in peaks:
            #     if (peak - 2) * WINDOWS < pos < (peak+2) * WINDOWS:
            #         inside = True
            # print p, inside
            # missed = missed or (not inside)


    

    # chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]

    # for offset, ratchet in zip(offsets.split(','), rs.split(',')):
    #     if len(offset) > 3 and offset[:4] == 'grav':
    #         if len(offset) > 4:
    #             count = int(offset[4:])
    #         else:
    #             count = 0
    #         if strand == '+':
    #             plt.axvline(int(ratchet) - int(start), linewidth=4, color='g')
    #         elif strand == '-':
    #             plt.axvline(int(end) - int(ratchet), linewidth=4, color='g')
    #     else:
    #         count = int(offset)

    #     if count:
    #         if strand == '+':
    #             plt.axvline(int(ratchet) - int(start), linewidth=2, color='r')
    #         elif strand == '-':
    #             plt.axvline(int(end) - int(ratchet), linewidth=2, color='r')

    if True:
        print chrom, start, end, strand, rs

        x = range(0, len(expression), WINDOWS)
        if len(mcmc) != len(x): x = range(0, len(expression)+1, WINDOWS)
        plt.plot(x, mcmc, linewidth=6, color= 'k')

        plt.plot(expression)

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
                plt.axvline(pos, ymax = .5, linewidth=4, color='r')

        # for gchrom, gstrand, grs in strad:
        #     if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
        #         if strand == '+':
        #             pos = grs - int(start)
        #         elif strand == '-':
        #             pos = int(end) - grs
        #         plt.axvline(pos, ymax = .33, linewidth=2, color='r')

        # for peak in peaks:
        #     peak.plot(len(expression))

        for i in range(30, len(seq) - 30):
            motif = seq[i - len(tp_pwm): i + len(fp_pwm)]

            score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)

            if score > .8:
 #               marker = '*' if increased(i, expression) else 'o'
                marker = 'o'
                plt.scatter([i], [max(expression)], marker = marker, linewidths = [(score - .8) * 100], c = color(score))

        plt.autoscale(tight = True)
        plt.show(block = False)
        a = raw_input("enter to continue")
        if a: plt.savefig('../data/coverage_example/' + a + '.png')
        plt.close()
print y, n
