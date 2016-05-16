from load_genome import load_genome, revcomp
from get_motifs import make_pwm, score_motif, get_min_score, get_max_score, search_for_motif
import sys
from math import exp

genome_seq = load_genome(open(sys.argv[4], 'r'))
fp_pwm, tp_pwm = make_pwm(sys.argv[3], genome_seq)

pwm = tp_pwm + fp_pwm
min_score, max_score = get_min_score(pwm), get_max_score(pwm)


class Intron:
        def __init__(self, chrom, five, strand):
		self.chrom = chrom
		self.five = five
		self.strand = strand

                # things that will be altered by EM
		self.per = []
		self.sjr = {}
		self.pi = None
		self.seq = None
                self.a_m = None
                self.scores = None
                self.counts = None
                self.assignments = None

                # fixed parameters
                self.a_s = 1.0
                self.mu = 1.0

	def key(self):
		return (self.chrom, self.five, self.strand)

	def add_sjr(self, three, count):
		i = abs(self.five - three)
		self.sjr[i] = count

	def add_per(self, three, five_insert):
		i = abs(self.five - three)
		self.per += [(i, five_insert)]
                
        def length(self):
                return abs(self.five - self._get_downstream_end())

	def _get_downstream_end(self):
		return max(map(lambda x: x[0], self.per) + self.sjr.keys())

        def _set_seq(self):
                three = self._get_downstream_end()
                start = min(self.five, three) - len(tp_pwm)
                end = max(self.five, three) + len(tp_pwm)
                self.seq = genome_seq[self.chrom][start:end]
                if self.strand == '-': self.seq = revcomp(self.seq)

        def _set_scores(self):
                self.scores = []
                for i in range(len(tp_pwm), len(self.seq) - len(fp_pwm)):
                        motif = self.seq[i-len(tp_pwm):i + len(fp_pwm)]
                        if motif[len(tp_pwm) - 2: len(tp_pwm) + 2] != 'AGGT': 
                                self.scores += [0]
                        else:
                                self.scores += [(score_motif(pwm, motif) - min_score) / (max_score - min_score)]

        def _set_a_m(self):
                assert self.scores != None
                best = max(self.scores)
                self.a_m = []
                for score in self.scores:
                        self.a_m += [self.a_s * max(0,  score - .8) * self.mu / best]

        def _set_pi(self):
                self.pi = [1 / float(len(self.scores))] * len(self.scores)


        def _e_step(self):
                """
                Here you update self.counts based on fractional read assignments.
                """
                self.counts = [0] * len(self.scores)

                # sjr are assigned deterministically
                for pos in self.sjr:
                        self.counts[pos] += self.sjr[pos]

                # per are assigned based on prob of generation
                for three, five_insert in self.per:
                        total = 0.0
                        for pos in range(max(0, three - 1200), three + 12):
                                if self.pi[pos] == 0: continue
                                inferred_insert = five_insert + three - pos
                                total += self._prob_insert(inferred_insert) * self.pi[pos]
                        if total == 0:
                                continue
                        for pos in range(max(0, three - 1200), three + 12):
                                if self.pi[pos] == 0: continue
                                inferred_insert = five_insert + three - pos
                                self.counts[pos] += (self._prob_insert(inferred_insert) * self.pi[pos]) / total

        def _assign_per(self):
                self.assignments = []
                for three, five_insert in self.per:
                        best = (-1, 0)
                        for pos in range(max(0, three - 1200), three + 12):
                                inferred_insert = five_insert + three - pos
                                p = self._prob_insert(inferred_insert) * self.pi[pos]
                                if p > best[1]:
                                        best = (pos, p)
                        if best[0] != -1:
                                self.assignments += [(three, five_insert, best[0], best[1])]

        def _m_step(self):
                """
                Here you update self.pi based on read assignments and prior
                
                do actual update in em so that you can compare old to new
                """
                pi = []
                for c, m in zip(self.counts, self.a_m):
                        pi += [max(0, c - self.a_s + m)]
                total = sum(pi)
                print total
                pi = map(lambda x: x / total, pi)
                return pi

        def _em(self):
                self._set_seq()
                self._set_scores()
                self._set_a_m()
                self._set_pi()

                done = False
                while not done:
                        self._e_step()
                        print sum(self.counts)
                        pi = self._m_step()
                        print sum(self.pi)
                        done = not any(abs(i-j) > .00001 for i, j in zip(pi, self.pi))
                        self.pi = pi
                
                self._assign_per()

        def _prob_insert(self, insert):
                """
                No normalization
                Keep numbers bigger...
                """
                u, o = 200.0, 30.0
                dev = insert - u
                if dev > 10 * o: return 0
                return exp(dev ** 2 / (-2*o**2))

        def i_to_g(self, pos):
                if self.strand == '+':
                        return self.five + end
                else:
                        return self.five - end

        def get_assignments(self):
                self._em()
                ends = {}
                for three, insert, prediction, score in self.assignments:
                        if prediction in ends:
                                ends[prediction] += 1
                        else:
                                ends[prediction] = 1
                out = []
                for end in ends:
                        out += [(self.chrom, self.five, self.i_to_g(three), self.strand, ends[end], self.pi[end], self.strand)]
                for i in self.sjr:
                        print self.i_to_g(three), self.sjr[i]
                return out


per_file = open(sys.argv[1], 'r')
introns = {}
for line in per_file:
	chrom, inner_left, inner_right, sample, five, strand = line.strip().split('\t')
	inner_left, inner_right, five = int(inner_left), int(inner_right), int(five)


	if strand ==  '+':
		five_insert = five - inner_left
		three = inner_right
	else:
		five_insert = inner_right - five
		three = inner_left

	key = (chrom, five, strand)

        if key not in introns:
                introns[key] = Intron(chrom, five, strand)
        introns[key].add_per(three, five_insert)



sjr_file = open(sys.argv[2])
for line in sjr_file:
	chrom, start, end, samples, count, strand = line.strip().split('\t')[:6]

	if strand == '+':
		five, three = int(start), int(end)
	else:
		five, three = int(end), int(start)

	key = (chrom, five, strand)

        if key in introns:
                introns[key].add_sjr(three, int(count))
                           
for intron in introns:
        assignments = introns[intron].get_assignments()
        for assign in assignments:
                print assign
                
                
