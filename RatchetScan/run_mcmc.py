from mcmc import mcmc
import sys

WINDOWS = 100
T = 5.0

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

data = open(sys.argv[1], 'r')
c = 0
with open(sys.argv[1]) as fp:
    for line in fp:
        start, end = line.split()[1:3]
        if int(end)-int(start) < WINDOWS*7: continue
        expression = [int(i) for i in line.strip().split('\t')[-1].split(",")]
        c += 1
        if line.split('\t')[-2] == "-":
            expression.reverse()

        samples = mcmc(expression, WINDOWS, T)

        print line.strip() + '\t' + ','.join(map(str, get_probs(samples, expression)))
