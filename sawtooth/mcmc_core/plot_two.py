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


data1 = open(sys.argv[1], 'r')
data2 = open(sys.argv[2], 'r')
graveley = open('../data/graveley.bed', 'r')

grav = []
for line in graveley:
    chrom, start, end, name, score, strand = line.strip().split('\t')
    if strand == '+':
        grav += [(chrom, strand, int(end))]
    else:
        grav += [(chrom, strand, int(start))]
c = 0
for line1, line2 in zip(data1, data2):
    c += 1
    print c
    chrom, start, end, offsets, rs, strand = line1.strip().split('\t')[:6]
    #if c < 65: continue   
    expression1 = [int(i) for i in line1.strip().split('\t')[-2].split(",")]
    expression2 = [int(i) for i in line2.strip().split('\t')[-2].split(",")]

    if line1.split('\t')[-3] == '-':
        expression1.reverse()
        expression2.reverse()

    mcmc1 = [float(i) for i in line1.strip().split('\t')[-1].split(",")]
    mcmc2 = [float(i) for i in line2.strip().split('\t')[-1].split(",")]

    for gchrom, gstrand, grs in grav:
        if gchrom == chrom and gstrand == strand and int(start) < grs < int(end):
            if strand == '+':
                plt.axvline(int(grs) - int(start), linewidth=4, color='g')
            elif strand == '-':
                plt.axvline(int(end) - int(grs), linewidth=4, color='g')


    
    # if line.split('\t')[5] == "-":
    #     expression.reverse()

    # print ' '.join(line.split('\t')[:6])

    if len(mcmc1) != len(range(0, len(expression1), WINDOWS)): continue

    plt.plot(range(0, len(expression1), WINDOWS), mcmc1, linewidth=6, color= 'k')
    plt.plot(range(0, len(expression2), WINDOWS), mcmc2, linewidth=4, color= 'r')



    # total = sum(expression)

    

    print chrom, start, end, rs, strand

    # # for offset, ratchet in zip(offsets.split(','), rs.split(',')):
    # #     if len(offset) > 3 and offset[:4] == 'grav':
    # #         if len(offset) > 4:
    # #             count = int(offset[4:])
    # #         else:
    # #             count = 0
    # #         if strand == '+':
    # #             plt.axvline(int(ratchet) - int(start), linewidth=4, color='g')
    # #         elif strand == '-':
    # #             plt.axvline(int(end) - int(ratchet), linewidth=4, color='g')
    # #     else:
    # #         count = int(offset)

    # #     if count:
    # #         if strand == '+':
    # #             plt.axvline(int(ratchet) - int(start), linewidth=2, color='r')
    # #         elif strand == '-':
    # #             plt.axvline(int(end) - int(ratchet), linewidth=2, color='r')

    plt.plot(expression1)
    plt.plot(expression2)
    plt.show(block = False)
    a = raw_input("enter to continue")
    plt.close()