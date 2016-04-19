from anno_jxns import Junctions
import pysam
from load_genome import load_genome, revcomp
from get_motifs import make_pwm, score_motif, get_min_score, get_max_score, search_for_motif
import sys

genome_seq = load_genome(open(sys.argv[3], 'r'))
fp_pwm, tp_pwm = make_pwm(sys.argv[2], genome_seq)

fp_min, fp_max = get_min_score(fp_pwm), get_max_score(fp_pwm)
tp_min, tp_max = get_min_score(tp_pwm), get_max_score(tp_pwm)

samfile = pysam.AlignmentFile(sys.argv[1])
jxns    = Junctions(sys.argv[2])

for read in samfile.fetch():
	if not read.is_read2: continue

	strand = (read.is_read1 == read.is_reverse)
	strand_str = '+' if strand else '-'
	chrom = samfile.getrname(read.reference_id)
	if chrom == 'Uextra': continue

	if strand and 1000 < read.pnext - read.pos:
		splices = jxns.search('chr' + chrom, strand_str, read.pos + 51, read.pnext)

		# check if consistent with any anno splice
		inner_left = read.pos + 51
		inner_right = read.pnext
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
				pwm = tp_pwm + fp_pwm
				start = read.pnext + (five - read.pos) - 300
				end = read.pnext + len(fp_pwm)
				seq = genome_seq[chrom][start:end]

				min_score = fp_min + tp_min
				max_score = fp_max + tp_max
				i, p = search_for_motif(pwm, seq, min_score, max_score)

				if p > .8:
					print chrom, five + 2, start + i + len(tp_pwm) + 1, strand_str, read.pos, read.pnext