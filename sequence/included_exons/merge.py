import sys
# python sequence/included_exons/merge.py ../data/five_included/motif.bed ../data/ten_included/motif.bed  ../data/twenty_included/motif.bed ../data/total_included/motif.bed

class Entry:
	def __init__(self, line):
		a = line.strip().split()
		self.chrom = a[0]
		self.start = int(a[1])
		self.end = int(a[2])
		self.position = a[3]
		self.body = int(a[4])
		self.strand = a[5]
		self.up = int(a[6])
		self.down = int(a[7])

five   = open(sys.argv[1])
ten    = open(sys.argv[2])
twenty = open(sys.argv[3])
total  = open(sys.argv[4])

for lines in zip(five, ten, twenty, total):
	entries = map(Entry, lines)


	bodies = map(lambda x: x.body, entries)
	ups    = ','.join(map(str, map(lambda x: x.up  , entries)))
	downs  = ','.join(map(str, map(lambda x: x.down, entries)))

	first = entries[0]

	bodies = ','.join(map(str, bodies))
	ups = ','.join(map(str, ups))
	downs = ','.join(map(str, downs))

	# bodies is montonically decreasing

	#if downs != '0,0,0,0': continue
	if ups == '0,0,0,0': continue

	print '\t'.join(map(str, [first.chrom, first.start, first.end, first.position, bodies, first.strand, ups, downs]))
