import sys
from intron_rs import IntronRS
from intron import Intron
import matplotlib.pyplot as plt

data = open(sys.argv[1], 'r')

rs = []
for line in data:
	intronrs = IntronRS(line)
	# if not intronrs.aggt(): continue
	# if not intronrs.motifs_above_thresh(.82): continue
	# if not intronrs.expression() > 10: continue
	rs += [Intron(intronrs)]

print len(rs)

introns = []
for r in rs:
	seen = False
	for intron in introns:
		if intron.compatible(r):
			intron.add(r)
			seen = True
			break
	if not seen:
		introns += [r]

print len(introns)

bed = open('../data/long_highly_expressed.bed')
lengths = []
for line in bed:
	start, end = line.strip().split('\t')[1:3]

	lengths += [int(end) - int(start)]


plt.title('Recursive Intron Lengths')
plt.hist(lengths, bins = 50)
plt.hist(map(lambda x: x.length(), introns), bins = 50, color = 'r')
plt.ylim([0, 100])
plt.savefig('../data/test/intron_lengths.png')
plt.show()
plt.title('Recursive Segment Lengths')
plt.hist(reduce(lambda x, y: x + y, map(lambda x: x.rs_lengths(), introns)), color = 'r', bins = 50)
plt.savefig('../data/test/rs_lengths.png')
plt.show()

bars = [0] * 10
for intron in introns:
	bars[intron.num_rs()-1] += 1
	if intron.num_rs() > 7:
		print intron.five(), intron.chrom
print bars

plt.title('Recursive Sites per Intron')
plt.bar(range(1, 11), bars, align = 'center')
plt.savefig('../data/test/num_rs.png')
plt.show()

plt.hist(map(lambda x: x.rs_length_dispersion(), introns))
plt.savefig('../data/test/dispersion.png')
plt.show()
