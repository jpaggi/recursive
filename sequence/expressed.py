import sys

for line in sys.stdin:
	chrom, start, end, name, sjr, strand = line.strip().split('\t')[:6]

	if int(sjr) > int(sys.argv[1]) and int(sys.argv[3]) > int(end) - int(start) >= int(sys.argv[2]):
			print line.strip()
