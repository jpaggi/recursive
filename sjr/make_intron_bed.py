from get_anno_splices import get_jxns

ss = get_jxns('../../data/anno.ss')
CONTEXT = 46

for site in ss:
	chrom, strand, five = site
	if strand == '+':
		three = max(ss[site])
		start, end = five + 1, three
	else:
		three = min(ss[site])
		start, end = three + 1, five

	if end - start < 500: continue

	start = start - CONTEXT
	end = end + CONTEXT

	print '\t'.join(map(str, [chrom[3:], start, end,'', 0, strand]))