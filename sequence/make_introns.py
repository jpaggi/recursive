from intron import Intron
import matplotlib.pyplot as plt


def main(entries):
	introns = []
	for line in open('../data/long_highly_expressed.bed'):
		chrom, start, end, name, counts, strand = line.strip().split('\t')[:6]
		introns += [Intron(chrom, int(start), int(end), strand, int(counts))]

	for entry in entries:
		chrom, rs, strand = entry.chrom, entry.rs, entry.strand
		for intron in introns:
			if intron.compatible(chrom, strand, rs):
				intron.add(rs)

	rs_introns = filter(lambda x: x.rs, introns)
	for intron in rs_introns:
		print intron

	return

	not_rs = filter(lambda x: not x.rs, introns)

	print "Introns over 400000 with no RS"
	for intron in not_rs:
		if intron.length() > 40000:
			print intron.chrom, intron.start, intron.endte

	

	plt.title('Recursive Intron Lengths')
	plt.hist(map(lambda x: x.length(), introns),    bins = 50)
	plt.hist(map(lambda x: x.length(), rs_introns), bins = 50, color = 'r')
	plt.ylim([0, 30])
	#plt.savefig('../data/test/intron_lengths.png')
	plt.show()
	plt.title('Recursive Segment Lengths')
	plt.hist(reduce(lambda x, y: x + y, map(lambda x: x.rs_lengths(), rs_introns)), color = 'r', bins = 50)
	#plt.savefig('../data/test/rs_lengths.png')
	plt.show()

	bars = [0] * 10
	for intron in rs_introns:
		bars[intron.num_rs()-1] += 1
	print bars

	plt.title('Recursive Sites per Intron')
	plt.bar(range(1, 11), bars, align = 'center')
	#plt.savefig('../data/test/num_rs.png')
	plt.show()

	plt.hist(map(lambda x: x.rs_length_dispersion(), rs_introns))
	#plt.savefig('../data/test/dispersion.png')
	plt.show()

if __name__ == '__main__':
	import sys
	from standard_table_reader import Entry
	main([Entry(line) for line in open(sys.argv[1], 'r')])
