import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
# import random

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