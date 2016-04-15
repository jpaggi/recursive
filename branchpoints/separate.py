def separate(ID, seq1, seq2, qual1, qual2):
	"""
	Input reads must both be forward strand (matching RNA)
	"""
	# read2 is broken
	# 1
	entry1 = "{}-{}\n".format(ID, seq2)
	entry2 = "{}-{}\n".format(ID, seq2)
	#2
	entry1 += seq2[-15:] + '\n'
	entry2 += seq1 + '\n'
	#3
	entry1 += '+\n'
	entry2 += '+\n'
	#4
	entry1 += qual2[-15:] + '\n'
	entry2 += qual1 + '\n'

	#read1 is broken
	#1
	entry1 += "{}-{}\n".format(ID, seq1)
	entry2 += "{}-{}\n".format(ID, seq1)
	#2
	entry1 += seq2 + '\n'
	entry2 += seq1[-15:] + '\n'
	#3
	entry1 += '+\n'
	entry2 += '+\n'
	#4
	entry1 += qual2 + '\n'
	entry2 += qual1[-15:] + '\n'

	return entry1, entry2
