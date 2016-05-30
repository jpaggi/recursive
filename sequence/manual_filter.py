import sys
import matplotlib.pyplot as plt
from standard_table_reader import Entry

expression = open('../data/all_merged.bed', 'r')
e = []
good = open(sys.argv[2], 'w')

for line in expression:
	chrom, start, end, name, count, strand, expression = line.strip().split('\t')
	start, end, count = int(start), int(end), int(count)
	expression = map(int, expression.split(','))
	e += [(chrom, start, end, strand, expression)]
print 'finsihed reading'


data = open(sys.argv[1], 'r')
for line in data:
	entry = Entry(line)

	if entry.manual or not entry.good_motif(.87):
		good.write(str(entry) + '\n')
		continue

	coverage = None
	for chrom, start, end, strand, expression in e:
		if chrom  == entry.chrom and start < entry.rs < end and strand == entry.strand:
			coverage = expression
			pos = entry.rs - start if strand == '+' else end - entry.rs
			break

	if coverage != None:
		print "{} {}".format(entry.chrom, entry.rs)
		print entry.motif_score
		if entry.strand == '-':
			coverage = coverage[::-1]
		plt.plot(coverage)
		plt.axvline(pos, c = 'r', linewidth=2)
		plt.show(block = False)

		a = raw_input('any key for yes')

		if a: entry.add_manual()
		good.write(str(entry) + '\n')

		plt.close()

	else:
		print 'no expression???'

