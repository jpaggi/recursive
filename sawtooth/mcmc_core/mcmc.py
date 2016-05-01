import numpy as np
from scipy import stats
from random import random, randrange, choice, gauss
from math import exp, log, ceil

class Fit:
    def __init__(self, x):
        self.x = x
        self.slopes = np.empty([len(x), len(x)], np.dtype('f16'))
        self.intercepts = np.empty([len(x), len(x)], np.dtype('f16'))
        self.devs = np.empty([len(x), len(x)], np.dtype('f16'))
        for i in range(len(x)):
            for j in range(i + 2, len(x)):
                self.slopes[i, j] = 0
                self.intercepts[i, j] = 0
                self.devs[i, j] = -1

    def set(self, i, j):
        if self.devs[i, j] == -1:
            slope, intercept, dev = weighted_regression(self.x, i, j)
            self.slopes[i, j] = slope
            self.intercepts[i, j] = intercept
            self.devs[i, j] = self.get_log_dev(i, j, slope, intercept)

    def get_log_dev(self, i, j, slope, intercept):
        avg = (slope*(1 + j - i)) / float(2)  + intercept
        if avg <= 0: return float('inf')
        dev = - sum(self.x[i:j]) * log((j - i) * avg)
        for p in range(i, j):
            height = slope*(p - i) + intercept
            if height <= 0: return float('inf')
            dev += self.x[p] * log(height)
        return dev


    def get_dev(self, i, j):
        self.set(i, j)
        return self.devs[i, j]

    def get_slope(self, i, j):
        self.set(i, j)
        return self.slopes[i, j]

    def get_int(self, i, j):
        self.set(i, j)
        return self.intercepts[i, j]

    def get_end_height(self, i, j):
        return self.get_int(i, j) + self.get_slope(i, j) * (j -i)

def weighted_regression(x, start, end):
    y = np.array(range(0, end - start))
    z = np.array(x[start:end])
    slope, intercept = np.polyfit(y, z, 1)
    if slope > 0:
        slope = 0
        intercept = sum(z) / float(len(z))
        res = sum((e - intercept) ** 2 for e in z) / float(intercept + .01)
        return slope, intercept, res
    i = 0
    while i < 5:
        w = [1 / float(max(1, intercept + slope * j)) for j in y]
        old_slope = slope
        slope, intercept = np.polyfit(y, z, 1, w = w)
        i += 1
        if abs(slope - old_slope) < .01: break

    w = [1 / float(max(1, intercept + slope * j)) for j in y]
    fit, residuals = np.polyfit(y, z, 1, w = w, full = True)[:2]
    res = residuals[0] if end - start > 2 else 0
    return fit[0], fit[1], res


def lin_regression(x, start, end):
    """
    calculate slope and y-intercept of data x
    beginning at start and ending at end
    
    requires that interval is at least length 2

    start is inclusive, end is exculsive
    """
    assert(end > start + 1)
    y = np.array(range(0, end - start))
    slope, intercept, r_value, p_value, std_err = stats.linregress(y, x[start:end]) # O(n) list slice
    return slope, intercept

def descending_regression(x, start, end):
    """
    Calculate slope and intercept of data x
    beginning at start and ending at end, 
    ensuring that the slope is non-positive

    requires that interval is at least length 2 

    start is inclusive, end is exclusive
    """
    #do normal linear regression and check if slope is non-positive
    slope, intercept = lin_regression(x, start, end)
    if slope <= 0:
        return slope, intercept
    # otherwise return slope = 0 and intercept = mean
    total = sum([x[i] for i in range(start, end)])
    return 0, total / float(end - start)

def calculate_squared_dev(x, start, end):
    """
    calculate the squared deviation of data x
    from the best fit line.
    No normalization not an r^2 value
    Currently uses descending regression

    start is inclusive, end is exclusive
    """
    slope, intercept = descending_regression(x, start, end)
    dev = 0
    for i in range(start, end):
        expect = intercept + slope * (i - start)
        dev += abs(x[i] - expect) ** 2
    return slope, intercept, dev

def score(state, fit, num):
    params = 2 * len(state) + 2
    if params == 2:
        rss = fit.get_dev(0, num - 1)
    else:
        rss = fit.get_dev(0, state[0])
        end_height = fit.get_end_height(0, state[0])
        for i in range(len(state) - 1):
            rss += fit.get_dev(state[i], state[i+1])
            if fit.get_int(state[i], state[i+1]) < end_height * 1.5:
                return float('inf')
            end_height = fit.get_end_height(state[i], state[i+1])
        rss += fit.get_dev(state[-1], num-1)
        if fit.get_int(state[-1], num-1) < end_height * 1.5:
            return float('inf')
    if rss == 0 or num == 0 or rss / float(num) <= 0:
        return - float('inf')
    return num * rss + 2 * params * log(num)

class State:
    def __init__(self, length, state = ()):
        assert type(state) == tuple
        self.length = length
        self.state = state

    def __eq__(self, other):
        return self.state == other.state

    def remove(self):
        if not self.state: return State(self.length, self.state)
        leaving = choice(self.state)
        new_state = tuple(x for x in self.state if x != leaving)
        return State(self.length, new_state)

    def add(self):
        entering = randrange(2, self.length-3)
        if any([(entering+j in self.state) for j in range(-1, 2)]):
            new_state = self.state
        else:
            new_state = tuple(sorted(self.state + (entering, )))
        return State(self.length, new_state)

    def move(self):
        if not self.state: return State(self.length, self.state)
        change = int(ceil(abs(gauss(0, 2))) * choice((-1, 1)))
        leaving = choice(self.state)
        entering = leaving + change

        if not 2 <= entering <= self.length - 3: 
            return State(self.length, self.state)

        temp_state = tuple(x for x in self.state if x != leaving)

        if any([(entering+j in temp_state) for j in range(-1, 2)]):
            new_state = self.state
        else:
            new_state = tuple(sorted(temp_state + (entering, )))
        return State(self.length, new_state)

    def get_state(self):
        return self.state


def mcmc(x, window, T):
    x = reduce_array(x, window)
    state = State(len(x))
    fit = Fit(x)
    old_score = score(state.get_state(), fit, len(x))
    samples = []
    for i in xrange(10000):
        transition_factor = 1.0
        cutoff = random()
        if cutoff < .4 and state:
            new_state = state.remove()
        elif cutoff < .6 and state:
            new_state = state.move()
        else:
            new_state = state.add()
            if len(new_state.get_state()) == 1: transition_factor = 2.0

        if state == new_state: continue

        new_score = score(new_state.get_state(), fit, len(x))

#        print state.get_state(), new_state.get_state(), old_score, new_score

        if old_score > new_score or new_score == - float('inf'):
            old_score = new_score
            state = new_state
        else:
            thresh = (exp(float(old_score - new_score) / T) / transition_factor)
            if thresh > random():
                old_score = new_score
                state = new_state

        if i > 1000 and not i % 50:
            samples += [state.get_state()]

    return samples

def reduce_array(array, windows):
    """
    create a new array with dimensionality reduced by a 
    factor of 'windows'

    currently cuts off partial window at the end
    
    averages entries in array up to this point.
    """
    length = int(len(array) / windows)
    out = np.zeros(length)
    for i in range(length):
        start = i * windows
        out[i] = sum([array[start + j] for j in range(windows)]) / windows 
    return out
