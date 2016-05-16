
class Kallisto:
	def __init__(self, line):
		self.ID, self.length, self.effective_length, self.counts, self.tpm = line.strip().split()

	def get_ID(self):
		return self.ID
		
	def coverage(self):
		return self.counts / self.effective_length

	def merge(self, other):
		self.counts += other.counts
		assert self.ID == other.ID

	def __str__(self):
		return self.ID + '\t' + self.coverage()

def get_all():
	genes = {}

	samples = ["{}_rep{}/abundances.tsv".format(time, rep) for time in [5, 10, 20, 'total'] for rep in ['A', 'B', 'C']]
	samples.remove('total_repC/abundances.tsv')

	for sample in samples:
		with open(sample, 'r') as f:
			for line in f:
				entry = Kallisto(line)

				if entry.get_ID() in genes:
					genes[entry.get_ID].merge(entry)
				else:
					genes[entry.get_ID] = entry
	return genes

if __name__ == '__main__':
	for entry in get_all():
		print entry
