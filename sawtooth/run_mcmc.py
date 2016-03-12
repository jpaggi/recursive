from mcmc import mcmc
import sys


data = open(sys.argv[1], 'r')
c = 0
for line in data:
    c += 1
    #if c < 10: continue
    expression = [int(i) for i in line.strip().split('\t')[-1].split(",")]
#    print line.strip().split('\t')[:-1]
    
    if line.split('\t')[-2] == "-": 
        expression.reverse()
    print mcmc(expression, 50)