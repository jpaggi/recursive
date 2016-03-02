from get_anno_splices import get_jxns
import sys

ss = get_jxns('../../data/anno.ss')

last_five, last_three, last_chrom, ratchets, offsets, last_strand = 0, 0, '', [], [], ''

total = 0

for line in sys.stdin:
	chrom, start, end, sample, offset, strand, seq1, seq2 = line.strip().split('\t')
	if not (seq1[-2:] == 'AG' and seq2[:2] == 'GT' and int(offset) > 0): continue

	if strand == '+':
		five, three = start, end
	else:
		five, three = end, start

	if last_chrom and last_chrom == chrom and last_strand == strand and last_five == five:
		ratchets += [three]
		offsets += [offset]
	else:
		if last_chrom:
			if last_strand == '+':
				last_start, last_end = last_five, str(max(ss[('chr' + last_chrom, last_strand, int(last_five) - 1)]))
			else:
				last_start, last_end = str(min(ss[('chr' + last_chrom, last_strand, int(last_five))])), last_five

			total += len(ratchets)
			print '\t'.join([last_chrom, last_start, last_end, ','.join(offsets), ','.join(ratchets), last_strand])


		last_chrom, last_five, last_three, offsets, last_strand = chrom, five, three, [offset], strand
		ratchets = [three]

if last_strand == '+':
	last_start, last_end = last_five, str(max(ss[('chr' + last_chrom, last_strand, int(last_five) - 1)]))
else:
	last_start, last_end = str(min(ss[('chr' + last_chrom, last_strand, int(last_five))])), last_five

print '\t'.join([last_chrom, last_start, last_end, ','.join(offsets), ','.join(ratchets), last_strand])

#print total + len(ratchets)
