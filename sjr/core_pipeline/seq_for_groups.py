import sys
from subprocess import Popen, PIPE
from load_genome import load_genome, revcomp
from get_motifs import make_pwm, score_motif, get_min_score, get_max_score

SEQ_LENGTH = 30

introns = open(sys.argv[1], 'r')
genome_seq = load_genome(open(sys.argv[2], 'r'))
out = open(sys.argv[3], 'w')

fp_pwm, tp_pwm = make_pwm(sys.argv[4], genome_seq)

fp_min, fp_max = get_min_score(fp_pwm), get_max_score(fp_pwm)
tp_min, tp_max = get_min_score(tp_pwm), get_max_score(tp_pwm)

for line in introns:
	chrom, start, end, sample, offsets, strand = line.strip().split('\t')[:6]
	if strand == '+':
		five, three = int(start), int(end)
	else:
		three, five = int(start), int(end)

	seq = genome_seq[chrom][three - SEQ_LENGTH:three + SEQ_LENGTH]

	if strand == '-': seq = revcomp(seq)

	fps = seq[SEQ_LENGTH:]
	tps = seq[:SEQ_LENGTH]

	fps_bits = (score_motif(fp_pwm, fps[:8]) - fp_min) / (fp_max - fp_min) 
	tps_bits = (score_motif(tp_pwm, tps[-15:]) - tp_min) / (tp_max - tp_min)

	if fps_bits > 0.8: print fps
	out.write('\t'.join([chrom, start, end, sample, offsets, strand, tps, fps, str(fps_bits), str(tps_bits)]) + '\n')
