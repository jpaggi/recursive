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
    return dev

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


def init_squared_dev(x, k):
    # n+1 x k arrays
    r = np.empty([len(x)+1, k], np.dtype('f8'))  # total squared deviation
    lookup = np.empty([len(x)+1, k], np.dtype('uint16')) # most recent subproblem used
    end = np.empty([len(x)+1, k], np.dtype('f8'))

    #initialize scores to infinity
    for i in range(len(x) + 1):
        for j in range(k):
            r[i, j] = float('inf')
    # dev always zero for zero bases
    r[0] = [0] * k
    lookup[0] = [0] * k

    #dev always zero for one base
    r[1] = [0] * k
    lookup[1] = [0] * k
    return r, lookup, end


def broken_regression(x, k, window):
    """
    use dynamic programming to compute max-likelihood segmentation for given region
    """
    # reduce dimension of array
    x = reduce_array(x, window)
    # n+1 x k arrays
    r, lookup, ends = init_squared_dev(x, k)

    # find deviations for no ratchet splicing up to all points
    for i in range(2, len(x) + 1):
        ends[i, 0] = end_height(x, 0, i)[1]
        lookup[i, 0] = 0
        r[i, 0] = 0

    # calculate remaining deviations recursively
    # i is ending index
    for i in range(2, len(x)+1): #start at 2 because zero for all at 1
        # j is beginning index of last spliced portion
        for j in range(i - 1):
            intercept, end = end_height(x, j, i)
            dev = calculate_squared_dev(x, j, i)

            # q is number of ratchet splices being considered
    
            for q in range(1, k):
                score = - intercept / max(1, ends[j][q-1]) + r[j][q - 1]
                if ends[j][q-1] >= intercept: score = float('inf')
                if score < r[i][q]:
                    r[i][q] = score
                    lookup[i][q] = j
                    ends[i][q] = end
    splices = back_track(r, lookup, k, window)
    return get_params(x, splices, window)

def back_track(r, l, k, window):
    """
    return list of length k with entries that are tuples:
    ( dev ^ 2 for k ratchets, [ratchet locations])
    """
    out = []
    for q in range(k):
        i = int(q)
        n = len(l) - 1
        entry = (r[n, i], [])
        while i >= 0 and n > 0:
            n = l[n, i]
            entry[1].append(n * window)
            i -= 1
        entry[1].reverse()
        out.append(entry)
    return out

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

def get_params(x, splices, window):
    out = []
    for splice in splices:
        params = []
        score = 0.0
        for i in range(len(splice[1])):
            start = splice[1][i] / window
            if i == len(splice[1]) - 1:
                end = len(x)
            else:
                end = splice[1][i+1] / window
            if end - start == 1:
                slope, intercept = 0, x[start]
            elif end - start > 1:
                slope, intercept = descending_regression(x, start, end)
                score += calculate_squared_dev(x, start, end)
            else:
                assert False
            params += [(slope / float(window), intercept)]
        #assert(abs(splice[0] - score) < .01)
        out.append((splice[0], splice[1], params))
    return out
