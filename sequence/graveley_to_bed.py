graveley = open('../data/recursiveintrons.bed')

for line in graveley:
	chrom, start, end, name, rs, strand  = line.strip().split('\t')

	for r in map(int, rs.split(',')):
		if strand == '+':
			begin, stop = start, r-1
		else:
			begin, stop = r, end
		print '\t'.join(map(str, [chrom, begin, stop, name, '.', strand]))