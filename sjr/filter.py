import sys
import matplotlib.pyplot as plt
from intron_rs import IntronRS

lengths = [0]
lengths_good_motif = [0]
lengths_great_motif = [0]
n, b, t = 0, 0, 0
for line in sys.stdin:
	intron = IntronRS(line)
	if not intron.expressed(): continue
	if not intron.aggt(): continue
	t += 1
	lengths += [intron.length()]
	if intron.length() >= 30000:
		pass
	if intron.recursive_index() > .05:
		print intron.motif_str()
		if intron.good_motifs():
			n += 1
		else:
			b += 1
	if intron.good_motifs():
		lengths_good_motif += [intron.length()]
	if intron.great_motifs():
		lengths_great_motif += [intron.length()]

print n, b, t
plt.title('Recursively Spliced Intron Lengths')
plt.xlabel('Intron Length')
plt.ylabel('Number of Recursive Sites')
plt.hist(lengths, bins = 50)
plt.hist(lengths_good_motif, color = 'r', bins = 50)
plt.hist(lengths_great_motif, color = 'g', bins = 50)

print len(lengths_good_motif)
print len(lengths_great_motif)

plt.ylim([0, 130])
plt.show()
