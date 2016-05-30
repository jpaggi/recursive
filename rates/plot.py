import matplotlib.pyplot as plt
import sys
from math import log

fxn = lambda x: log(x) if sys.argv[2] == 'log' else x

class MISO:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom    = a[0]
		self.start    = int(a[1])
		self.end      = int(a[2])
		self.strand   = a[3]
		self.rs_start = int(a[4])
		self.rs_end   = int(a[5])
		self.num      = int(a[6])
		self.total_rs = int(a[7])
		self.mean_psi = map(fxn, map(lambda x: float(x) if float(x) > 0 else 0.01, a[8].split(',')))

	def five(self):
		return self.mean_psi[:3]

	def ten(self):
		return self.mean_psi[3:6]

	def twenty(self):
		return self.mean_psi[6:9]

	def rs_length(self):
		return self.rs_end - self.rs_start

	def total(self):
		return self.mean_psi[9:]

f = open(sys.argv[1])
f.readline()
f.readline()

lengths = []
five_min_psi = []
for line in f:
	miso = MISO(line)

	plt.scatter([fxn(5)]*3, miso.five())
	plt.scatter([fxn(10)]*3, miso.ten())
	plt.scatter([fxn(20)]*3, miso.twenty())
	plt.scatter([fxn(60)]*2, miso.total())

	print "length: {}, pos {} of {}".format(miso.rs_length(), miso.num, miso.total_rs)
	plt.show(block = False)
	a = raw_input('enter')
	plt.close()


	# lengths +=  [miso.rs_length()]
	# five_min_psi += [sum(miso.five()) / 3.0]

# plt.scatter(lengths, five_min_psi)
# plt.show()
