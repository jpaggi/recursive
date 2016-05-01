from anno_jxns import Junctions
import pysam
#from load_genome import load_genome, revcomp
#from get_motifs import make_pwm, score_motif, get_min_score, get_max_score, search_for_motif
import sys

#genome_seq = load_genome(open(sys.argv[3], 'r'))
#fp_pwm, tp_pwm = make_pwm(sys.argv[2], genome_seq)

#pwm = tp_pwm + fp_pwm
#min_score, max_score = get_min_score(pwm), get_max_score(pwm)

samfile = pysam.AlignmentFile(sys.argv[1])
jxns    = Junctions(sys.argv[2])

for read in samfile.fetch():
	if not read.is_read2: continue
	if read.reference_id != read.next_reference_id: continue

	strand = (read.is_read1 == read.is_reverse)
	strand_str = '+' if strand else '-'
	chrom = samfile.getrname(read.reference_id)
	if chrom == 'Uextra': continue
	if chrom == 'dmel_mitochondrion_genome': continue

	if strand and 1000 < read.pnext - read.pos:
		inner_left = read.pos + 51
		inner_right = read.pnext
		splices = jxns.search('chr' + chrom, strand_str, inner_left, inner_right)

		# check if consistent with any anno splice
		consistent = False
		for five, three in splices:
			cassette = five >= inner_left and three <= inner_right
			if cassette and five - inner_left + inner_right - three < 450:
				consistent = True
		if consistent: continue

		fives = {}
		for five, three in splices:
			if five in fives:
				fives[five] = max(three, fives[five])
			else:
				fives[five] = three


		for five in fives:
			three = fives[five]
			if five - read.pos < 300 and read.pnext < three:
				# start = read.pnext + (five - read.pos) - 300
				# end = read.pnext + len(fp_pwm)
				# seq = genome_seq[chrom][start:end]
				# i, p = search_for_motif(pwm, seq, min_score, max_score)

				# if p > .8:
				# 	print chrom, five + 2, start + i + len(tp_pwm) + 1, strand_str, read.pos, read.pnext
				print '\t'.join(map(str, [chrom, inner_left, inner_right, 'PE', five + 1, strand_str]))

	elif not strand and 1000 < read.pos - read.pnext:
		inner_left  = read.pnext
		inner_right = read.pos
		splices = jxns.search('chr' + chrom, strand_str, inner_left, inner_right)

		consistent = False
		for five, three in splices:
			cassette = three >= inner_left and five <= inner_right
			if cassette and three - inner_left + inner_right - five < 450:
				consistent = True
		if consistent: continue

		fives = {}
		for five, three in splices:
			if five in fives:
				fives[five] = min(three, fives[five])
			else:
				fives[five] = three


		for five in fives:
			three = fives[five]
			if inner_right - five < 300 and inner_left > three:
				# start = inner_left - len(fp_pwm)
				# end = inner_left + 300 - inner_right + five
				# seq = revcomp(genome_seq[chrom][start:end])

				# i, p = search_for_motif(pwm, seq, min_score, max_score)

				# if p > .8:
				# 	print chrom, end - i - len(tp_pwm) - 2, five, strand_str, read.pnext, read.pos
				print '\t'.join(map(str, [chrom, inner_left, inner_right, 'PE', five, strand_str]))
