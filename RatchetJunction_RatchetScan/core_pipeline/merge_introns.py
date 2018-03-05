import sys
from intron_rs import IntronRS

cur = IntronRS(sys.stdin.readline())
for line in sys.stdin:
	rs = IntronRS(line)
	if cur.same(rs):
		cur.merge(rs)
	else:
		print cur
		cur = rs
print cur
