import matplotlib.pyplot as plt
from broken_regression import broken_regression
import sys



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



for line in sys.stdin:
    chrom, start, end, name, a, strand, expression, fit = line.strip().split("\t")
    start = int(start)
    end = int(end)
    expression = map(lambda x: int(x), expression.split(","))
    fit = eval(fit)
    if strand == "-": 
        expression.reverse()
    for splices in fit:
        plt.plot(plot_fit(splices[1], splices[2], len(expression)))
        plt.plot(expression)
        plt.show()
