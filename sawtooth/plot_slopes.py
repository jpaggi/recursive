import matplotlib.pyplot as plt
import sys

data = open(sys.argv[1], 'r')


import numpy as np
from scipy import stats

def lin_regression(x, start, end):
    """
    calculate slope and y-intercept of data x
    beginning at start and ending at end
    
    requires that interval is at least length 2

    start is inclusive, end is exculsive
    """
    assert(end > start + 1)
    y = np.array(range(0, end - start))
    slope, intercept, r_value, p_value, std_err = stats.linregress(y, x[start:end]) # O(n) list slice
    return slope, intercept

count = 0
for line in data:
    #print line
    expression = [int(i) for i in line.strip().split('\t')[-1].split(",")]

    total = sum(expression)

    chrom, start, end, offsets, rs, strand = line.strip().split('\t')[:6]

    if strand == '+':
        begin = int(start)
        for ratchet in sorted(map(int, rs.split(','))):
            print lin_regression(expression, begin - int(start), ratchet - int(start))
            begin = ratchet
        print lin_regression(expression, begin - int(start), int(end) - int(start))

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

    if total / float(int(end) - int(start)) > 1:
    	pass	
    print '\t'.join(line.split('\t')[:-1])
    
    if line.split('\t')[-2] == "-": 
        expression.reverse()
    plt.plot(expression)
    plt.show(block = False)
    a = raw_input("enter to continue")
    plt.close()