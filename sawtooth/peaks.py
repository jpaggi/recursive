WINDOWS = 100
from random import randrange
class Peak:
    def __init__(self, start, end, best_prob, genome_space = False):
        self.best_prob = best_prob
        if not genome_space:
            self.start = (start-1) * WINDOWS
            self.end = end * WINDOWS
        else:
            self.start, self.end = start, end

    def inside(self, pos):
        return self.start <= pos <= self.end

    def close(self, pos, dist):
        return self.start - dist <= pos <= self.end + dist

    def length(self):
        return self.end - self.start

    def random(self, length):
        if length - self.length() <= 0:
            return Peak(0, length, self.best_prob, True)
        start = randrange(length - self.length() - 1)
        end = start + self.length()
        return Peak(start, end, self.best_prob,  True)

    def plot(self, length):
        start_p = self.start / float(length)
        end_p = self.end / float(length)
        plt.axvline(self.start, c = 'm')
        plt.axvline(self.end, c = 'm')

    def mergable(self, end):
        end = end * WINDOWS
        return end < self.end + 1000

    def merge(self, end, best_prob):
        self.end = end * WINDOWS
        self.best_prob = max(best_prob, self.best_prob)

    def __str__(self):
        return "{},{}".format(self.start, self.end)

def get_peaks(probs, z, THRESH):
    """
    Given probs in reduced space
    Return Peaks in genome space

    Currently missing peaks that go all the way to end!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    """
    if THRESH == 0: return [Peak(0, len(probs), max(probs) / z)]
    peaks = []
    peak_start = None
    best_prob = 0
    for i in xrange(len(probs)):
        p = probs[i] / float(z)
        # make new peak
        if peak_start != None and p < THRESH:
            if peaks and peaks[-1].mergable(i):
                peaks[-1].merge(i, best_prob)
            else:
                peaks += [Peak(peak_start, i, best_prob)]
            peak_start = None
            best_prob = 0
        elif peak_start == None and p >= THRESH:
            peak_start = i
            best_prob = max(best_prob, p)
    return peaks
