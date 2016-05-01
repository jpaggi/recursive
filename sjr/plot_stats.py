import matplotlib.pyplot as plt
import sys

data = open(sys.argv[1], 'r')


ratio, five, three  = [], [], []
for line in data:
	chrom, start, end, sample, rs, strand, sjr, total , seq1, seq2, fps, tps = line.strip().split('\t')

	start, end, sjr, total, fps, tps = int(start), int(end), int(sjr), int(total), float(fps), float(tps)

	if sjr / float(total + 1) < .9:
		ratio += [sjr / float(total + 1)]
		five += [fps]
		three += [tps]


plt.scatter(ratio, three)
plt.show()