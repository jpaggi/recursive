import sys
from standard_table_reader import Entry
import matplotlib.pyplot as plt
from random import randrange
from math import log, exp, pi, sqrt

def FILTER(arr, thresh, s = .2):
	return (len(filter(lambda x: thresh-s <= x <= thresh+s, arr)) + 1) / float(len(arr) + 1)


def INCREASE(coverage, pos):
	LENGTH = 500
	up = coverage[max(0, pos - LENGTH): pos]
	down = coverage[pos: min(pos + LENGTH, len(coverage))]
	up_avg = sum(up) / float(len(up))
	dn_avg = sum(down) / float(len(down))

	increase = min(8, (dn_avg+1) / (up_avg+1))

	return increase if increase < 2 else float(int(increase))

def GAUSS(arr):
	true_mean = sum(arr) / float(len(arr))
	true_var  = sum(map(lambda x: (x - true_mean) ** 2, arr)) / float(len(arr))
	return lambda x: exp((- (x - true_mean) **2) / (2 * true_var)) / sqrt(2 * pi * true_var)


random_motifs = map(float, open('../data/aggt_motif_scores.csv').read()[:-1].split(','))
true_motifs = []

expression = open('../data/coverage/all_merged_masked.bed', 'r')
e = []
for line in expression:
	chrom, start, end, name, count, strand, expression = line.strip().split('\t')
	start, end, count = int(start), int(end), int(count)
	expression = map(int, expression.split(','))
	e += [(chrom, start, end, strand, expression)]
print 'finsihed reading'


#############################################################################################
########    Calculate increases in expression at observed graveley sites     ################
########        and get random increases 5 x upsampled                       ################
true_increase = []
random_increase = []
data = open(sys.argv[1], 'r')
for line in data:
	if line[0] == 'C': continue
	entry = Entry(line)
	if not (sum(entry.junc) or entry.saw_score): continue
	coverage = None
	for chrom, start, end, strand, expression in e:
		if chrom  == entry.chrom and start < entry.rs < end and strand == entry.strand:
			coverage = expression
			pos = entry.rs - start if strand == '+' else end - entry.rs
			r_positions = [randrange(500, len(expression) - 500) for i in range(5)]
			break

	if coverage != None:
		if entry.strand == '-': coverage = coverage[::-1]

		if entry.grav:
			true_increase += [INCREASE(coverage, pos)]
			true_motifs += [entry.motif_score]

			for r_pos in r_positions:
				random_increase += [INCREASE(coverage, r_pos)]

data.close()

# For increase distribution use empirical dist with some rounding
#  ... to give necessary smoothing.
# Fit Gaussian to motif distributions, needed for smoothing
t_m = GAUSS(true_motifs)
r_m = GAUSS(random_motifs)


plt.hist(true_increase, bins = 50)
plt.title('True Increase Distribution')
plt.show()
plt.hist(random_increase, bins = 50)
plt.title('Random Increase Distribution')
plt.show()

plt.hist(true_motifs, bins = 50)
plt.title('True Motif Distribution')
plt.show()
plt.hist(random_motifs, bins = 50)
plt.title('Random Motif Distribution')
plt.show()

#########################################################################################################

random_sawtooth = [float(line.strip().split()[4]) for line in open('../data/mcmc_peak_calls/all_merged_masked_random.bed')]

true_sawtooth = []
data = open(sys.argv[1], 'r')
for line in data:
	if line[0] == 'C': continue
	entry = Entry(line)
	if not (sum(entry.junc) or entry.saw_score): continue

	if entry.grav: true_sawtooth += [entry.saw_score]
data.close()

#########################################################################################################

print 'ROC Curves if using strict cutoff in any one parameter.'
x, y = [], []
for i in [0,.01, .02, .03, .04, .05, .06, .07, .08, .09, .1, .15, .2]:
	x += [len(filter(lambda x: x >= i, random_sawtooth)) / float(len(random_sawtooth))]
	y += [len(filter(lambda x: x >= i, true_sawtooth)) / float(len(true_sawtooth))]

plt.plot(x, y)
plt.title('Sawtooth ROC')
plt.show()

x, y = [], []
for i in [0, .1, .2, .3, .4, .5, .6, .7, .8, .82, .83, .84, .85, .86, .87, .88, .89, .9, .92, .95]:
	x += [len(filter(lambda x: x >= i, random_motifs)) / float(len(random_motifs))]
	y += [len(filter(lambda x: x >= i, true_motifs)) / float(len(true_motifs))]

plt.plot(x, y)
plt.title('Motif ROC')
plt.show()

x, y = [], []
for i in [0, .5, .75, 1, 1.25, 1.5, 2, 3, 4, 6, 10]:
	x += [len(filter(lambda x: x >= i, random_increase)) / float(len(random_increase))]
	y += [len(filter(lambda x: x >= i, true_increase)) / float(len(true_increase))]
plt.plot(x, y)
plt.title('Coverage Increase ROC')
plt.show()

############################################################################################################

print 'Get ROC curve combining expression increase and motif scores'
true_log = []
data = open(sys.argv[1], 'r')
for line in data:
	if line[0] == 'C': continue
	entry = Entry(line)
	if not (sum(entry.junc) or entry.saw_score): continue

	coverage = None
	for chrom, start, end, strand, expression in e:
		if chrom  == entry.chrom and start < entry.rs < end and strand == entry.strand:
			coverage = expression
			pos = entry.rs - start if strand == '+' else end - entry.rs
			r_pos = randrange(500, len(expression) - 500)
			break

	if coverage != None:
		if entry.strand == '-':
			coverage = coverage[::-1]

		if entry.grav:
			increase = INCREASE(coverage, pos)

			p_increase_plus  = FILTER(true_increase, increase)
			p_increase_minus = FILTER(random_increase, increase)

			p_motif_plus  = t_m(entry.motif_score)
			p_motif_minus = r_m(entry.motif_score)

			log_increase = log(p_increase_plus) - log(p_increase_minus)
			log_motif    = log(p_motif_plus) - log(p_motif_minus)

			log_total = log_motif + log_increase
			true_log += [log_total]

data.close()

random_log = []
for i, increase in enumerate(random_increase):
	motif = random_motifs[i]

	# p_increase_plus  = FILTER(true_increase, increase, .2)
	# p_increase_minus = FILTER(random_increase, increase, .2)

	# p_motif_plus  = FILTER(true_motifs, motif)
	# p_motif_minus = FILTER(random_motifs, motif)

	p_increase_plus  = FILTER(true_increase, increase)
	p_increase_minus = FILTER(random_increase, increase)

	p_motif_plus  = t_m(motif)
	p_motif_minus = r_m(motif)

	log_increase = log(p_increase_plus) - log(p_increase_minus)
	log_motif    = log(p_motif_plus) - log(p_motif_minus)

	log_total = log_motif +  log_increase

	random_log += [log_total]



x, y = [], []
cut = 10
for i in [-float('inf')]+range(-15, 15):
	x += [len(filter(lambda x: x > i, random_log)) / float(len(random_log))]
	y += [len(filter(lambda x: x > i, true_log)) / float(len(true_log))]

	print i, y[-1], x[-1]
	if y[-1] < .1 and i < cut:
		print 'cutoff', i, y[-1]

plt.plot(x, y)
plt.show()


#####################################################################################################################
# APPLY TO ALL DATA!!!!!!!!!!!

t_g, f_g, t_n, f_n = 0, 0, 0, 0
true_log = []

data = open(sys.argv[1], 'r')
out = open(sys.argv[2], 'w')
for line in data:
	if line[0] == 'C': continue
	entry = Entry(line)

	coverage = None
	for chrom, start, end, strand, expression in e:
		if chrom  == entry.chrom and start < entry.rs < end and strand == entry.strand:
			coverage = expression
			pos = entry.rs - start if strand == '+' else end - entry.rs
			r_pos = randrange(500, len(expression) - 500)
			break

	if coverage != None:
		if entry.strand == '-':
			coverage = coverage[::-1]

		increase = INCREASE(coverage, pos)


		p_increase_plus  = FILTER(true_increase, increase)
		p_increase_minus = FILTER(random_increase, increase)

		p_motif_plus  = t_m(entry.motif_score)
		p_motif_minus = r_m(entry.motif_score)

		log_increase = log(p_increase_plus) - log(p_increase_minus)
		log_motif    = log(p_motif_plus) - log(p_motif_minus)

		log_total =  log_motif + log_increase

		true_log += [log_total]

		if log_total > -1 and entry.grav:
			t_g += 1
		elif log_total > -1 and not entry.grav:
			t_n += 1
		elif entry.grav:
			f_g += 1
		else:
			f_n += 1

		entry.log_or = log_total
		out.write(str(entry) + '\n')
data.close()

print 'Perfomance using cutoff of -1:', t_g, f_g, t_n, f_n

plt.hist(true_log, bins = 50)
plt.show()




