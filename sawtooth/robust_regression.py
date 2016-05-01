from sklearn import linear_model, datasets
from matplotlib import pyplot as plt
import numpy as np
import sys
WINDOWS = 1




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


data = open(sys.argv[1], 'r')
c = 0
for line in data:
    c += 1
    if c < 1: continue
    expression = [int(i) for i in line.strip().split('\t')[6].split(",")]
    if line.split('\t')[-2] == "-":
        expression.reverse()

    y = reduce_array(expression, WINDOWS)
    X = np.array([range(len(y))]).transpose()

    print X

    

    model = linear_model.LinearRegression()

    
    model.fit(X, y)

    # Robustly fit linear model with RANSAC algorithm
    model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression())
    try:
        model_ransac.fit(X, y)
    except ValueError:
        continue
    inlier_mask = model_ransac.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)

    # Predict data of estimated models
    line_X = np.array(range(len(y)))
    line_y = model.predict(line_X[:, np.newaxis])
    line_y_ransac = model_ransac.predict(line_X[:, np.newaxis])

    # Compare estimated coefficients
    print("Estimated coefficients (true, normal, RANSAC):")
    print(model.coef_, model_ransac.estimator_.coef_)

    plt.plot(X[inlier_mask], y[inlier_mask], '.g', label='Inliers')
    plt.plot(X[outlier_mask], y[outlier_mask], '.r', label='Outliers')
    plt.plot(line_X, line_y, '-k', label='Linear regressor')
    plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
    plt.show()
