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

x = []
y = []

events = []

for line in data:
    chrom, start, end, name, a, strand, expression, fit = line.strip().split('\t')
    fit = eval(fit)
    expression = map(lambda x: int(x), expression.split(","))
    if strand == "-":
        expression = expression[::-1]
    start = int(start)
    end = int(end)
    length = end - start
    if length < 2000: continue
    cpb = sum(expression) / length
    for i in [5]:#range(len(fit)):
        if float(fit[0][0]) == 0: continue
        x += [i]
        y += [float(fit[i][0]) / float(fit[0][0])]
        events += [(y[-1], expression, fit)]
events = sorted(events, key = lambda x: x[0])
for event in events:
    score, expression, fit = event
    splices = fit[5]
    
    plt.plot(plot_fit(splices[1], splices[2], len(expression)))
    plt.plot(expression)
    plt.show(block = False)
    a = raw_input("enter to continue")
    plt.close()
