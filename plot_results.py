import matplotlib.pyplot as plt
from broken_regression import broken_regression


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



data = open("exons_cov.bed", 'r')
c = 0
for line in data:
    expression = [int(i) for i in line.strip().split('\t')[-1].split(",")]
#    print line.strip().split('\t')[:-1]
    
    if line.split('\t')[-2] == "-": 
        expression.reverse()
    if 5000 > int(line.split('\t')[2]) - int(line.split('\t')[1]) > 500:
        try:
            fit = broken_regression(expression, 5, 20)
        except ValueError:
            # error when long run with same expression
            # these are generally strange looking so I think
            # it is okay to just discount them
            continue
        print line.strip() + '\t' + str(fit)
#        print c
        c += 1
#        for splices in fit:
#            print splices
#            plt.plot(plot_fit(splices[1], splices[2], len(expression)))
#            plt.plot(expression)
#            plt.show()
