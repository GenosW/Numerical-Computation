import numpy as np

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

if __name__ == "__main__":
    x = np.array([0,1,2])
    y = np.array([0,2,4])
    d = divdiff(x, y, x.size)
    print(d)
    ytest = horner(d, x, 10)
    print(ytest)
