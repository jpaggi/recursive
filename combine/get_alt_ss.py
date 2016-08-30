

data = open('../data/rs_jxns/all.bed', 'r')

sites = {} # RS to 5'SS

for line in data:
	chrom, start, end, rs, counts, strand, spanning = line.strip().split()

	five = int(start) if strand == '+' else int(end)
	counts = sum(map(int, counts.split(',')))

	RS_STR = "{}:{}:{}".format(chrom, rs, strand)
	FIVE_STR = "{}:{}:{}".format(chrom, five, strand)

	entry = (FIVE_STR, counts, int(spanning))

	if RS_STR in sites:
		sites[RS_STR].append(entry)
	else:
		sites[RS_STR] = [entry]

for site in sites:
	entries = sites[site]
	if len(entries) > 1:
		entries = ';'.join(["{},{},{}".format(*entry) for entry in entries])
		print "{}\t{}".format(site, entries)