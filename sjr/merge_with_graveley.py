import sys

introns = {}


def my_index(query, template):
	if query in template:
		return template.index(query)
	else:
		return -1

for line in open(sys.argv[1], 'r'):
	chrom, start, end, igv, ratchets, strand = line.strip().split('\t')
	ratchets = map(int, ratchets.split(','))

	if strand == '+':
		five, three = int(start), int(end)
	else:
		five, three = int(end), int(start)

	introns[(chrom, five, strand)] = (three, ['grav'] * len(ratchets), ratchets)

count = 0
for line in open(sys.argv[2], 'r'):
	chrom, start, end, offsets, ratchets, strand = line.strip().split('\t')
	ratchets = map(int, ratchets.split(','))
	if strand == '+':
		five, three = int(start) + 1, int(end)
		ratchets = map(lambda x: x + 1, ratchets)
	else:
		five, three = int(end), int(start)

	if (chrom, five, strand) in introns:
		for ratchet, offset in zip(ratchets, offsets.split(',')):
			pos = my_index(ratchet, introns[(chrom, five, strand)][2])
			if pos != -1:
				count += 1
				introns[(chrom, five, strand)][1][pos] = 'grav' + offset
			else:
				introns[(chrom, five, strand)][2].append(ratchet)
				introns[(chrom, five, strand)][1].append(offset)
	else:
		introns[(chrom, five, strand)] = (three, offsets.split(','), ratchets)


for intron in introns:
	chrom, five, strand = intron
	three, offsets, ratchets = introns[intron]

	if strand == '+':
		start, end = five, three
	else:
		start, end = three, five

	print '\t'.join(map(str, [chrom, start, end, ','.join(map(str, offsets)), ','.join(map(str, ratchets)), strand]))
