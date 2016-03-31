import numpy as np
from scipy import stats
from random import random, randrange, choice, gauss
from math import exp, log

def calculate_squared_dev(x, start, end, slope, intercept):
    """
    calculate the squared deviation of data x
    from the best fit line.
    No normalization not an r^2 value
    Currently uses descending regression

    start is inclusive, end is exclusive
    """
    dev = 0
    for i in range(start, end):
        expect = intercept + slope * (i - start)
        dev += (x[i] - expect) ** 2
    return dev

def calculate_rss(x, state, slope, intercepts):
    bounds = [0] + state + [len(x)]
    rss = 0
    for begin, end, intercept in zip(bounds[:-1], bounds[1:], intercepts):
        local = calculate_squared_dev(x, begin, end, slope, intercept)
        #print x, slope, intercepts, state
        #assert local <= calculate_squared_dev(x, begin, end, slope, intercept - 5)
        assert local <= calculate_squared_dev(x, begin, end, slope, intercept + 1)
        rss += local
    return rss

def remove(l, r):
    out = []
    for el in l:
        if el != r:
            out += [el]
    return out

def calculate_slope(state, x, NUM_MEMO):
    num = 0
    bounds = [0] + state + [len(x)]
    for begin, end in zip(bounds[:-1], bounds[1:]):
        if (begin, end) in NUM_MEMO:
            num += NUM_MEMO[(begin, end)]
        else:
            contrib = 0
            s_avg = (end + begin + 1) / float(2)
            x_avg = sum(x[begin:end]) / float(end - begin)
            for i in xrange(begin, end):
                contrib += (i - s_avg) * (x_avg - x[i])

            NUM_MEMO[(begin, end)] = contrib
            num += contrib

    denom = 0
    for begin, end in zip(bounds[:-1], bounds[1:]):
        s_avg = (end - begin + 1) / float(2)
        denom += (s_avg * (s_avg + 1) * (2 * s_avg + 1)) / float(3)

    return num / denom

def calculate_intercepts(state, x, m):
    bounds = [0] + state + [len(x)]
    y = []
    for begin, end in zip(bounds[:-1], bounds[1:]):
        s_avg = (end - begin + 1) / float(2)
        x_avg = sum(x[begin: end]) / float(end - begin)

        y += [x_avg - m * s_avg]
    return y

def score(state, x, NUM_MEMO):
    params = len(state) + 2
    num = len(x)
    
    slope = calculate_slope(state, x, NUM_MEMO)
    slope = min(0, slope)

    intercepts = calculate_intercepts(state, x, slope)

    # print x
    # print state, slope, intercepts

    bounds = [0] + state + [len(x)]
    for int1, int2, begin, end in zip(intercepts[:-1], intercepts[1:], bounds[:-1], bounds[1:]):
        if int1 + slope * (end - begin + 2) > int2:
            return float('inf')

    rss = calculate_rss(x, state, slope, intercepts)

    if rss == 0 or num == 0 or rss / float(num) <= 0:
        return - float('inf')
    return num * log(rss / float(num)) + 2 * params * log(num)

def mcmc(x, window):

    x = reduce_array(x, window)

    state = []

    NUM_MEMO = {}

    SCORE_MEMO = {}
    old_score = score(state, x, NUM_MEMO)

    samples = []
    for i in xrange(10000):
        SCORE_MEMO[tuple(state)] = old_score
        cutoff = random()
        if cutoff < .4 and state:
            new_state = remove(state, choice(state))
        elif cutoff < .6 and state:
            change = choice(state)
            new = int(gauss(change, 5))
            if (new not in state) and 0 < new < len(x):
                new_state = sorted(remove(state, change) + [new])
        else:
            add = randrange(1, len(x) - 1)
            if add in state: continue
            new_state = sorted(state + [add])

        new_score = SCORE_MEMO[tuple(new_state)] if tuple(new_state) in SCORE_MEMO else score(new_state, x, NUM_MEMO)

        if old_score > new_score or (new_score == - float('inf') or (old_score != - float('inf') and exp(old_score - new_score) > random())):
            old_score = new_score
            state = new_state

        if i > 1000 and not i % 10:
            samples += [state]

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
