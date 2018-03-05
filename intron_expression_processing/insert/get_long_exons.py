import sys

for line in sys.stdin:
	if line[0] == '#': continue
	if line[0] == '>': break
	event, start, end = line.strip().split()[2:5]

	if event != 'exon': continue
	if int(end) - int(start) < int(sys.argv[1]): continue

	print line.strip()
