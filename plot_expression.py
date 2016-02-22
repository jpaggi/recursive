import matplotlib.pyplot as plt
import sys

data = open(sys.argv[1], 'r')
for line in data:
    print line
    expression = [int(i) for i in line.strip().split('\t')[-1].split(",")]
    
    if line.split('\t')[-2] == "-": 
        expression.reverse()

    plt.plot(expression)
    plt.show(block = False)
    a = raw_input("enter to continue")
    plt.close()
