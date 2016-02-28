from random import randrange
from matplotlib import pyplot

class RNA:
    gene = [(0, 500), (1500, 2000)] # characterized by CDSs

    def __init__(self, extension):
        self.pos = extension
        self.degrade = 0

    def extend(self):
        if self.pos < self.gene[-1][1]:
            self.pos += 1
        else:
            self.degrade += 1
    def get_coverage(self):
        if self.pos <= 1000:
            return [(0, self.pos)]
        elif self.degrade >  100000:
            return []
        elif self.pos <= 1500:
            return [(0, 500), (1000, self.pos)]
        else:
            return [(0, 500), (1500, self.pos)]

def plot_coverage(rna):
    coverage = [0] * 2000
    for r in rna:
        for region in r.get_coverage():
            for i in range(region[0], region[1]):
                coverage[i] += 1
    maximum = float(max(coverage))
    return map(lambda x: x / maximum, coverage)

def plot_ideal_rna_seq(rna):
    coverage = [0] * 2000
    read_length = 50
    for r in rna:
        for region in r.get_coverage():
            for i in range(region[0], region[1] - read_length):
                for j in range(i, i + read_length):
                    coverage[j] += 1
    maximum = float(max(coverage))
    return map(lambda x: x / maximum, coverage)

def plot_ideal_pe_rna_seq(rna):
    coverage = [0] * 2000
    read_length = 50
    insert_length = 200
    for r in rna:
        for region in r.get_coverage():
            for i in range(region[0], region[1] - insert_length):
                for j in range(i, i + read_length) + range(i + insert_length - read_length, i + insert_length):
                    coverage[j] += 1
    maximum = float(max(coverage))
    return map(lambda x: x / maximum, coverage)

#initialize molecules
rna = []
for i in range(0):
    rna += [RNA(randrange(0, 2000))]
for i in range(3000):
    for r in rna:
        r.extend()
    rna += [RNA(0)] * 1



pyplot.plot(plot_coverage(rna))

pyplot.plot(plot_ideal_rna_seq(rna))

pyplot.plot(plot_ideal_pe_rna_seq(rna))


pyplot.show()
