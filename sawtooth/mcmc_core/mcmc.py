import numpy as np
from scipy import stats
from random import random, randrange, choice, gauss
from math import exp, log

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
        dev += abs(x[i] - expect) ** 22
    return slope, intercept, dev

def end_height(x, start, end):
    """
    calculate the squared deviation of data x
    from the best fit line.
    No normalization not an r^2 value
    Currently uses descending regression

    start is inclusive, end is exclusive
    """
    slope, intercept = descending_regression(x, start, end)
    return intercept, intercept + slope * (end - start)


def init(x):
    slopes = np.empty([len(x), len(x)], np.dtype('f8'))
    intercepts = np.empty([len(x), len(x)], np.dtype('f8'))
    devs = np.empty([len(x), len(x)], np.dtype('f8'))
    for i in range(len(x)):
        for j in range(i + 2, len(x)):
            slope, intercept, dev = calculate_squared_dev(x, i, j)
            slopes[i, j] = slope
            intercepts[i, j] = intercept
            devs[i, j] = dev
    return slopes, intercepts, devs

def remove(l, r):
    out = []
    for el in l:
        if el != r:
            out += [el]
    return out


def score(state, slopes, intercepts, dev):
    params = 2 * len(state) + 2
    num = len(intercepts[0])
    if params == 2:
        rss = dev[0, num - 1]
    else:
        rss = dev[0, state[0]]
        for i in range(len(state) - 1):
            rss += dev[state[i], state[i+1]]
        rss += dev[state[-1], num-1]


    if rss == 0 or num == 0 or rss / float(num) <= 0:
        return - float('inf')
    return num * log(rss) + 2 * params * log(num)

def mcmc(x, window):

    x = reduce_array(x, window)

    state = []
    slopes, intercepts, devs = init(x)

    old_score = score(state, slopes, intercepts, devs)

    samples = []
    for i in xrange(1000000):
        cutoff = random()
        if cutoff < .4 and state:
            new_state = remove(state, choice(state))
        elif cutoff < .6 and state:
            change = choice(state)
            new = int(gauss(change, 5))
            if (new not in state) and 0 <= new < len(x):
                new_state = sorted(remove(state, change) + [new])
        else:
            add = randrange(len(x - 1))
            if add in state: continue
            new_state = sorted(state + [add])

        new_score = score(new_state, slopes, intercepts, devs)
        #print x
        #print new_state, new_score, old_score

        thresh = exp(float(old_score - new_score))

        if old_score > new_score or (new_score == - float('inf') or (old_score != - float('inf') and thresh > random())):
            old_score = new_score
            state = new_state

        if i > 10000 and not i % 10:
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
