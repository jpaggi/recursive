"""
package to calculate optimal "broken" linear regressions

primary function called should be broken_regression(data, k, window)

where data is an array of integers, 
 k is a non-negative integer representing the number of splices made,
 and window is the size of bins used in averaging / size reduction


returns a length k array of tuples

each tuple corresponds to a (squared_dev, [list of splice sites]) pair

where the ith entry uses i ratchet sites indexed by zero-indexed first base of the window
where the break is made. Basically a list of 5'ss
"""

import numpy as np
from scipy import stats
from random import random, randrange, choice
from math import exp

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
        dev += (x[i] - expect) ** 2
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
    return slopes, intercepts, dev

def remove(l, r):
    out = []
    for el in l:
        if el != r:
            out += [el]
    return out


def score(state, slopes, intercepts, dev):
    return 1

def mcmc(x, window):

    x = reduce_array(x, window)

    state = []
    slopes, intercepts, devs = init(x)

    old_score = score(state, slopes, intercepts, devs)

    samples = []

    for i in xrange(100):
        if random() < .5 and state:
            new_state = remove(state, choice(state))
        else:
            add = randrange(len(x - 1))
            if add in state: continue
            new_state = sorted(state + [add])

        new_score = score(new_state, slopes, intercepts, devs)

        if exp(old_score - new_score) < random():
            old_score = new_score
            state = new_state

        if not i % 10:
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
