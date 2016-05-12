#META
# take a file of sjr in individual intron-centric format
# and an output directory
# plot following figures and save to output directory

# Goal is to do this in a completely automated way
# should be able to press play once and get all statistics
# might make sense to leave each as a script though?

import sys
from intron_rs import IntronRS

data = open(sys.argv[1], 'r')
directory = sys.argv[2]
rs = [IntronRS(line) for line in data]

print len(rs)


#overlaps with graveley study
# ... Number of graveley vs. number of samples
# numbers by number of samples
#  ... y axis is number of sjr
#  ... x axis is number of samples
# spectrum of motif strength???
import graveley_across_replicates
grav, novel = graveley_across_replicates.main(rs, True, directory)
print grav, novel

# numbers by length of intron
#  ... y axis is number of sjr
#  ... x axis is length of introns
# spectrum of motif strength
import intron_lengths
lengths, good_motif_lengths, great_motif_lengths = intron_lengths.main(rs, True, directory)
#print lengths, good_motif_lengths, great_motif_lengths

# numbers by recursive index
#  ... y axis is CDF of number of sjr
#  ... x axis is cutoff recursive index
# spectrum of motif strength... lol
import recursive_index
recursive_index.main(rs, True, directory)

# numbers by motif strength
# ... y axis is CDF of number of sjr
# ... x axis is cutoff of motif strangth
# spectrum of number of samples
