from mcmc import mcmc
import sys
import matplotlib.pyplot as plt
WINDOWS = 100

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
for line in data:
    c += 1
    #if c < 197: continue   
    expression = [int(i) for i in line.strip().split('\t')[-1].split(",")]

    if len(expression) < 500: continue
    
    if line.split('\t')[-2] == "-":
        expression.reverse()

    # samples = mcmc(expression, WINDOWS)

    # print line.strip() + '\t' + ','.join(map(str, get_probs(samples, expression)))

    # # for ratchet in samples[-1]:
    # #     print WINDOWS * ratchet
    # #     plt.axvline(WINDOWS * ratchet, linewidth=2, color='m')

    # #plt.plot(range(0, len(expression), WINDOWS), get_probs(samples, expression), linewidth=6, color= 'k')


    # total = sum(expression)

    # chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]

    # print rs

    # # for offset, ratchet in zip(offsets.split(','), rs.split(',')):
    #     print offset, ratchet
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

    dist = [0] * (max(expression) + 1)
    for e in expression:
        dist[e] += 1

    plt.plot(dist)
    plt.show(block = False)
    a = raw_input("enter to continue")
    plt.close()




