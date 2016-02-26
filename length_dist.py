import matplotlib.pyplot as plt
import sys

number = [0] * 100000

for line in open(sys.argv[1], 'r'):
	chrom, start, end, five, three, strand = line.strip().split('\t')
	length = abs(int(five) - int(three))
	if length < 100000: 
		number[length] += 1
	else:
		print length


plt.plot(number)
plt.show()