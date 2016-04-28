import sys
import matplotlib.pyplot as plt
WINDOWS = 100

sjr = open(sys.argv[1], 'r')
mcmc_file = open(sys.argv[1], 'r')

def get_prob(pos, start, mcmc, z):
    window = abs(pos - start) / WINDOWS
    score = 0
    for i in range(-5, 10):
        try:
            score = max(score, mcmc[window + i])
        except IndexError: pass
    return score / float(z)

# read mcmc into memory
for line in mcmc_file:
    chrom, start, end, name, sjr, strand, expression, samples = line.strip().split('\t')
    expression = map(int, expression.split(","))
    samples = map(float, samples.split(","))
    start, end, sjr = int(start), int(end), int(sjr)
    z = max(expression)

    if strand == "-":
        expression.reverse()

    if len(samples) != len(range(0, len(expression), WINDOWS)): continue
    begin = int(start) if strand == '+' else int(end)
    missed = False

    for r in map(int, rs.split(',')):
        p = get_prob(r, begin, samples, z)
        if p > .01:
            y += 1
        else:
            n += 1
            missed = True
    if missed:
        plt.plot(range(0, len(expression), WINDOWS), mcmc, linewidth=6, color= 'k')
        total = sum(expression)
        chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]

        for offset, ratchet in zip(offsets.split(','), rs.split(',')):
            if len(offset) > 3 and offset[:4] == 'grav':
                if len(offset) > 4:
                    count = int(offset[4:])
                else:
                    count = 0
                if strand == '+':
                    plt.axvline(int(ratchet) - int(start), linewidth=4, color='g')
                elif strand == '-':
                    plt.axvline(int(end) - int(ratchet), linewidth=4, color='g')
            else:
                count = int(offset)

            if count:
                if strand == '+':
                    plt.axvline(int(ratchet) - int(start), linewidth=2, color='r')
                elif strand == '-':
                    plt.axvline(int(end) - int(ratchet), linewidth=2, color='r')
        plt.plot(expression)
        plt.show(block = False)
        print line
        print "{}:{}-{}".format(chrom, start, end), rs
        a = raw_input("enter to continue")
        if a:
            problem.write(line)
        plt.close()

print y, n