import matplotlib.pyplot as plt
import sys

class Peak:
	def __init__(self, line):
		a = line.strip().split('\t')
		self.chrom  = a[0]
		self.start  = int(a[1])
		self.end    = int(a[2])
		self.prob   = float(a[3])
		self.score  = float(a[4])
		self.strand = a[5]


peak_scores = []
for line in open(sys.argv[1]):
	peak = Peak(line)
	peak_scores += [peak.score]

r_peak_scores = []
for line in open(sys.argv[2]):
	peak = Peak(line)
	r_peak_scores += [peak.score]

print len(filter(lambda x: x > .06, peak_scores))
print len(filter(lambda x: x > .06, r_peak_scores))

plt.hist(peak_scores, bins = 100)
plt.hist(r_peak_scores, bins = 100, color = 'r')
plt.show()