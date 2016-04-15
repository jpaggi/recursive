import sys
import matplotlib.pyplot as plt
WINDOWS = 100

data = open(sys.argv[1], 'r')
if len(sys.argv) > 2:
    problem = open(sys.argv[2], 'w')

y, n = 0, 0
def get_prob(pos, start, mcmc, z):
    window = abs(pos - start) / WINDOWS
    score = 0
    for i in range(-5, 10):
        try:
            score = max(score, mcmc[window + i])
        except IndexError: pass
    return score / float(z)

for line in data: 
    expression = [int(i) for i in line.strip().split('\t')[-2].split(",")]
    z = max(expression)
    mcmc = [float(i) for i in line.strip().split('\t')[-1].split(",")]

    chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]

    if line.split('\t')[5] == "-":
        expression.reverse()

    if len(mcmc) != len(range(0, len(expression), WINDOWS)): continue

    begin = int(start) if strand == '+' else int(end)

    missed = False
    for r in map(int, rs.split(',')):
        p = get_prob(r, begin, mcmc, z)

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

