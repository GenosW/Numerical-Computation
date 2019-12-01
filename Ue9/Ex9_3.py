import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
# import random

def newton1D(x, f, df):
    return x - f(x)/df(x)

def iterateNewton1D(x0, f, df, maxIteration, eps, bound=1000, debug=0):
    x = newton1D(x0,f,df)
    xList = [x0,x]
    iterations = 0
    for i in list(range(0,maxIteration)):
        if debug: print("Iteration: ", i+1,"\nx= ", x)
        x = newton1D(x,f,df)
        xList.append(x)
        if (abs(f(x)) < eps):
            iterations = i+1
            break
        if (abs(x-x0) >= bound):
            print("Out of bound! Terminating Newton-Raphson method.")
            iterations = i+1
            break
    return xList, iterations

def f(x):
    assert(isinstance(x,np.matrix))
    ret = np.matrix([[(x[0]-x[1])**2],[0]], dtype=float)
    return ret

if __name__ == "__main__":
    A = np.matrix([[2,1],[1,2]], dtype=float)
    b = np.matrix('0;3', dtype=float)
    x = np.matrix('0;0', dtype=float)
    xList = []
    eps = 1e-2
    x = np.linalg.solve(A, b)
    for i in range(0,10):
        iterations = i+1
        x = np.linalg.solve(A, b+eps*f(x))
        xList.append(np.array(x))
    print(iterations)
    print(xList)
    # print(A)
    # print(b)
    # ret = f(x)
    # print(ret)