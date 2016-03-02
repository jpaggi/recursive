import matplotlib.pyplot as plt
import sys

data = open(sys.argv[1], 'r')


introns = []
if False:#len(sys.argv) > 2:
	interest = open(sys.argv[2], 'r')

	for line in interest:
		chrom, start, end = line.strip().split('\t')[:3]
		introns.append((chrom, start, end))

novel = []
excluded = []
if len(sys.argv) > 3:
	putative = open(sys.argv[3], 'r')

	for line in putative:
		chrom, start, end, a, offsets, strand, seq1, seq2 = line.strip().split('\t')
		if seq1[-2:] == 'AG' and seq2[:2] == 'GT' and int(offsets) > 1:
			novel.append((chrom, start, end))
		else:
			excluded.append((chrom, start, end))


count = 0
for line in data:
    #print line
    expression = [int(i) for i in line.strip().split('\t')[-1].split(",")]

    total = sum(expression)

    chrom, start, end, igv, rs, strand = line.strip().split('\t')[:6]

    if ((chrom, start, end) not in introns) and introns: continue

    for i in rs.split(','):
    	if strand == '+':
    		plt.axvline(int(i) - int(start), linewidth=4, color='r')
    	else:
    		plt.axvline(int(end) - int(i), linewidth=4, color='r')

    for event in novel:
    	if chrom == event[0]:
    		if strand == '+' and int(start) < int(event[2]) < int(end):
    			plt.axvline(int(event[2]) - int(start), linewidth=2, color='g')
    		elif strand == '-' and int(start) < int(event[1]) < int(end):
    			plt.axvline(int(end) - int(event[1]), linewidth=2, color='g')
    for event in excluded:
    	if chrom == event[0]:
    		if strand == '+' and int(start) < int(event[2]) < int(end):
    			plt.axvline(int(event[2]) - int(start), linewidth=2, color='m')
    		elif strand == '-' and int(start) < int(event[1]) < int(end):
    			plt.axvline(int(end) - int(event[1]), linewidth=2, color='m')

    if total / float(int(end) - int(start)) > 1:
    	pass	
    print '\t'.join(line.split('\t')[:-1])
    
    if line.split('\t')[-2] == "-": 
        expression.reverse()
    plt.plot(expression)
    plt.show(block = False)
    a = raw_input("enter to continue")
    plt.close()