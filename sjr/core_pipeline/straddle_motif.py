from load_genome import load_genome, revcomp
from get_motifs import make_pwm, score_motif, get_min_score, get_max_score, search_for_motif
import sys

genome_seq = load_genome(open(sys.argv[3], 'r'))
fp_pwm, tp_pwm = make_pwm(sys.argv[2], genome_seq)

pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)

data = open(sys.argv[1], 'r')


for line in data:
	chrom, inner_left, inner_right, sample, five, strand = line.strip().split('\t')
	inner_left, inner_right, five = int(inner_left), int(inner_right), int(five)

	if strand ==  '+':
		five_insert = five - inner_left
		assert -5 <= five_insert <= 300

		start = inner_right + five_insert - 300
		end = inner_right + 5 + len(fp_pwm)
		seq = genome_seq[chrom][start:end]
		i, p = search_for_motif(pwm, seq, min_score, max_score, AG = True)
		rs = start + i + len(tp_pwm)
		if i == -1: continue
		print '\t'.join(map(str, [chrom, five, rs, sample, str(p), strand, ','.join(map(str, [inner_left, inner_right]))]))
	else:
		five_insert = inner_right - five
		assert -5 <= five_insert <= 300

		start = inner_left - len(fp_pwm)
		end = inner_left + 300 - five_insert
		seq = revcomp(genome_seq[chrom][start:end])
		i, p = search_for_motif(pwm, seq, min_score, max_score, AG = True)
		rs = end - i - len(tp_pwm)
		if i == -1: continue
		print '\t'.join(map(str, [chrom, rs, five, sample, p, strand, ','.join(map(str, [inner_left, inner_right]))]))
