graveley = open('../data/recursiveintrons.bed')

for line in graveley:
	chrom, start, end, name, rs, strand  = line.strip().split('\t')
	start = int(start)

	for r in map(int, rs.split(',')):
		if strand == '+':
			begin, stop = start-1, r-1
		else:
			begin, stop = r, end
		if chrom == '2R' and strand == '+' and stop == 13193490: stop += 1 
		print '\t'.join(map(str, [chrom, begin, stop, name, '.', strand]))