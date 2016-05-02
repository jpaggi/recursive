import sys
from get_motifs import *
from load_genome import *
WINDOWS = 100


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

straddle = open('../data/all_straddle_seq.bed', 'r')

strad = []
for line in straddle:
    chrom, start, end, name, score, strand, seq1, seq2 = line.strip().split('\t')

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

probs = []

for line in data:

    chrom, start, end, offsets, rs, strand, expression, mcmc = line.strip().split('\t')
    expression = map(int, expression.split(","))
    mcmc = map(float, mcmc.split(','))
    start, end = int(start), int(end)
    
    seq = genome_seq[chrom][int(start):int(end)]
    if strand == "-":
        expression.reverse()
        seq = revcomp(seq)

    for gchrom, gstrand, grs in grav:
        if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
            if strand == '+':
                pos = grs - int(start)
            elif strand == '-':
                pos = int(end) - grs
            probs += [get_prob(pos, mcmc, max(expression))]

    for gchrom, gstrand, grs in nov:
        if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
            if strand == '+':
                pos = grs - int(start)
            elif strand == '-':
                pos = int(end) - grs

    for gchrom, gstrand, grs in strad:
        if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
            if strand == '+':
                pos = grs - int(start)
            elif strand == '-':
                pos = int(end) - grs

    for i in range(30, len(seq) - 30):
        motif = seq[i - len(tp_pwm): i + len(fp_pwm)]

        score = (score_motif(pwm, motif) - min_score) / (max_score - min_score)

print probs

