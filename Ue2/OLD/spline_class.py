import numpy as np
import math as m
import matplotlib.pyplot as plt

def partition(a,b,N = 10):
    h = (b-a)/N
    t = list(map(lambda j: a+j*h,range(N))) # creates partition knots, range(N) gives j=0...N-1
    return t, h

def polKnots(t1, t2, n):
    x = list(map(lambda i: t1+i*(t2-t1)/n,range(n)))
    return x

def horner(d, knots, x):
    if isinstance(d, list):
        n = len(d)
    elif isinstance(d, np.ndarray):
        n = d.size - 1
    y = d[n]
    for j in range(n, -1, -1):
        y = d[j] + (x - knots[j])*y
    return y

def divdiff(x, y, n):
    f = np.zeros((n,n), dtype = float)
    f[:,0] = y
    for j in range(1,n):
        for i in range(n-j-1):
            f[i,j] = (f[i+1,j-1] - f[i,j-1])/(x[i+j]-x[i])
    return f[0,:]

def f(x):
    return m.sin(x*m.pi)

class spline:
    'My custom spline class :)'
    count = 0
   
    def __init__(self, x, y, N=10, p=3, r=2, name = 'My Spline'):
        self.x = x
        self.y = y
        self.N = N
        self.p = p
        self.r = r
        self.name = name
        spline.count += 1

    def __del__(self):
        class_name = self.__class__.__name__
        print(class_name, "destroyed")
    

if __name__ == "__main__":
    # x = np.array([0,1,2])
    # y = np.array([0,2,4])
    n = 2
    t, h = partition(0,4)
    knots = []
    fvals = []
    coeffs = []
    for (t1, t2) in zip(t[0:-1], t[1:]):
        polPoints = polKnots(t1,t2,n)
        knots.extend(polPoints)
        y = list(map(f, polPoints))
        fvals.extend(y)
        d = list(divdiff(polPoints, y, n))
        coeffs.extend(d)
    knots
    # print(knots, fvals)
    # print(coeffs)
    fig, ax = plt.subplots()
    ax.plot(knots, fvals)
    plt.show()

    # print(d)
    # ytest = horner(d, x, 10)
    # print(ytest)
