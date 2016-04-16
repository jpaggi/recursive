import sys
import matplotlib.pyplot as plt

data = open(sys.argv[1], 'r')
c = 0
for line in data:
    expression = [int(i) for i in line.strip().split('\t')[-1].split(",")]
#    print line.strip().split('\t')[:-1]
    
    if line.split('\t')[-2] == "-": 
        expression.reverse()
    slopes = []
    for i in range(len(expression) - 1):
        slopes += [expression[i+1] - expression[i]]

    slopes = sorted(slopes)

    avg_slope = sum(slopes) / float(len(slopes))

    slopes[len(slopes) / 4: (len(slopes) * 3) / 4]

    middle = sum(slopes[len(slopes) / 4: (len(slopes) * 3) / 4]) / float(len(slopes) / 2)

    print avg_slope, middle

