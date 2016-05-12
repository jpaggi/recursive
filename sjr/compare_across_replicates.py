import sys

for line in open(sys.argv[1], 'r'):
	chrom, start, end, reps, offsets, strand, fps, tps = line.strip().split('\t')

	reps = reps.split(',')
	reps = map(lambda x: x.split('_rep'), reps)
	offsets = int(offsets)

	five = 0
	ten = 0
	twenty = 0
	for rep in reps:
		if rep[0] == '5':
			five += 1
		elif rep[0] == '10':
			ten += 1
		else:
			twenty += 1

	if abs(int(start) - int(end)) < 5000: continue

	if ten + five + twenty == 3: print line.strip()
