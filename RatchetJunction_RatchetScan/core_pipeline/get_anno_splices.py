
def get_jxns(anno_file):
	"""
	Return a dictionary of splicing events present in the anno_file

	The anno file should have exons named well 'exon' and can be in
	either gff or gtf format

	The dictionary is (chrom, strand, five) -> [<associated 3'ss>, ...]
	"""
	
	sites =  open(anno_file, 'r')

	out = {}
	for line in sites:
		chrom, start, end, strand = line.strip().split('\t')[:7]
		if strand == '+':
			five, three = int(start), int(end)
		else:
			three, five = int(start), int(end)

		key = (chrom, strand, five)

		if key in out:
			out[key].append(three)
		else:
			out[key] = [three]
	return out
