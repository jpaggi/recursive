import sys
import matplotlib.pyplot as plt


def plot_fit(splices, params, length):
    y = []
    for i in range(len(splices)):
        start = splices[i] 
        if i == len(splices) - 1:
            end = length
        else:
            end = splices[i+1]
        slope, y_int = params[i]
        for j in range(end - start):
            y += [slope * j + y_int]
    return y

data = open(sys.argv[1], 'r')
SPLICES = int(sys.argv[2])
MIN_SIZE = int(sys.argv[3])

x = []
y = []

events = []

def f_test(rss1, rss2, p1, p2, n):
    num = (rss1 - rss2) / float((p2 - p1))
    denom = rss2 / float(n - p2 + 1)
    return num / denom
                               


for line in data:
    chrom, start, end, name, a, strand, expression, fit = line.strip().split('\t')
    fit = eval(fit)
    expression = map(lambda x: int(x), expression.split(","))
    if strand == "-":
        expression = expression[::-1]
    start = int(start)
    end = int(end)
    length = len(expression)
    if length < MIN_SIZE: continue
    cpb = sum(expression) / length
    if float(fit[0][0]) == 0: continue
    #ratio = float(fit[SPLICES][0]) / float(fit[0][0])
    rss1 = fit[0][0]
    rss2 = fit[SPLICES][0]
    if rss2 == 0: continue
    ratio = f_test(rss1, rss2, 2, SPLICES * 2, length)
    events += [(ratio, expression, fit, chrom, start, end)]

events = sorted(events, key = lambda x: x[0])
for event in events[::-1]:
    score, expression, fit, chrom, start, end = event
    splices = fit[SPLICES]
    print chrom, start, end
    plt.plot(plot_fit(splices[1], splices[2], len(expression)))
    plt.plot(expression)
    plt.show(block = False)
    a = raw_input("enter to continue")
    plt.close()
