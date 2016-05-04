import sys
import matplotlib.pyplot as plt
from repeats import Repeats
import numpy as np

def reduce_array(array, windows):
    """
    create a new array with dimensionality reduced by a 
    factor of 'windows'

    currently cuts off partial window at the end
    
    averages entries in array up to this point.
    """
    length = int(len(array) / windows)
    out = np.zeros(length)
    for i in range(length):
        start = i * windows
        out[i] = sum([array[start + j] for j in range(windows)]) / windows 
    return out

repeats = Repeats('../data/downloaded/repeats.txt')

data = open(sys.argv[1], 'r')
for line in data:
    chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]
    expression = [int(i) for i in line.strip().split('\t')[6].split(",")]
    start, end = int(start), int(end)
    if strand == "-":
        continue
        expression.reverse()

    # for b, s in repeats.get_repeats(chrom, start, end):
    #     if strand == '+':
    #         begin, stop = b - start, s - start
    #     else:
    #         begin, stop = end - s, end - b
    #     for pos in range(max(0, begin), min(stop, len(expression))):
    #         plt.axvline(pos, ymin = .9, c = 'm')

    expression = repeats.mask(chrom, start, end, strand, expression)

    #plt.plot(expression)

    #x = reduce_array(expression, 100)

    #plt.plot(range(0, len(expression), 100)[:-1], x)

    print '\t'.join([chrom, str(start), str(end), offsets, rs, strand, ','.join(map(str, expression))])

    #plt.show(block = False)
    #a = raw_input("enter to continue")
    #plt.close()
