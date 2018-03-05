import matplotlib.pyplot as plt
from math import exp

random = filter(lambda x: x < 50000, map(int, open('../data/peak_dists/random_dists.csv').read().strip().split(',')))
real = filter(lambda x: x < 50000, map(int, open('../data/peak_dists/real_dists.csv').read().strip().split(',')))
print 1+max(map(abs, random))


CDF = lambda dists: [len(filter(lambda x: abs(x) <= i, dists))/float(len(dists)) for i in range(1+max(map(abs, dists)))]


half = 500
k = -6 / float(half)

def logistic(x, k, x_0):
    if abs(x - x_0) > 2000: return 0
    return 1 / float(1 + exp(-k * (x - x_0)))


l = [logistic(i, k, 500) for i in range(1+max(map(abs, random)))]

random_cdf = CDF(random)
real_cdf = CDF(real)

PLOT = plt.plot # plot.semilogx

PLOT(random_cdf)
PLOT(real_cdf)
PLOT(l)

plt.savefig('../data/dists.png')
plt.show()