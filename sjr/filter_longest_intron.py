import sys

class Entry:
	def __init__(self, line):
		self.data = line.strip().split('\t')
		self.line = line
		self.length = int(self.data[1]) - int(self.data[2])

	def same(self, other):
		return self.data[0] == other.data[0] and self.data[3] == other.data[3]

	def longer(self, other):
		return self.length > other.length


entries = [Entry(sys.stdin.readline())]

for line in sys.stdin:
	entry = Entry(line)

	if not entry.same(entries[-1]):
		entries.append(entry)
	elif entry.longer(entries[-1]):
		entries[-1] = entry

for entry in entries:
	print entry.line.strip()

