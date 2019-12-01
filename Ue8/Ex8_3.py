import numpy as np
import matplotlib.pyplot as plt
import os
import random

def f(x):
    return np.sin(x)

def getRandKnots(a, b, N):
    knots = [0]*N
    for k in range(N):
        knots[k] = a+(b-a)*random.random()
    return np.array(knots)

class mySpecPolynom:
    def __init__(self, deg):
        self.a = [0]*deg
        self.deg = deg
        self.knots = []

    def setCoeffs(self, coeffs):
        self.a = coeffs

    def solveLS(self, xi, yi):
        self.knots = xi
        A = []
        b = []
        for row in range(0,2):
            A.append([np.sum(np.power(xi, i+row*2)) for i in (2,4)])
            b.append(np.dot(np.power(xi,1+row*2),yi))
        A = np.array(A)
        # print(A)
        # print(b)
        u = np.linalg.solve(A, b)
        # print(u)
        self.setCoeffs(u)

    def eval(self, x):
        return np.multiply(self.a[0],x) + np.multiply(self.a[1],np.power(x,3))


if __name__ == "__main__":
    scriptpath = os.path.dirname(__file__)
    pic = os.path.join(scriptpath, "Ex8_3_LS.png")
    Ns = [2**n for n in range(2,11)]
    # N = 2
    # Xi = getRandKnots(-1/N, 1/N, N)
    # Yi = np.sin(Xi)
    print("Test: ")
    # print(np.sum(np.power(Xi,2)))

    # foo = mySpecPolynom(2)
    # foo.solveLS(Xi, Yi)
    # print(foo.a)
    # print(foo.eval([-1,0,1]))
    print("Real stuff happens now")
    print("----------------------------------")
    Foos = []
    for N in Ns:
        foo = mySpecPolynom(N)
        Xi = getRandKnots(-1/N, 1/N, N)
        Yi = np.sin(Xi)
        foo.solveLS(Xi, Yi)
        Foos.append(foo)
        print(N,foo.a)

    reso = 10
    h = 1/reso
    x = np.arange(-1,1,h)
    # fig, ax = plt.subplots()
    for foo, sign in zip(Foos[-3:],('x','o','p')):
        plt.plot(x, foo.eval(x),sign, label=str(foo.deg))
    plt.plot(x, f(x), label="sin")
    # plt.plot(x, np.cos(x), label="cos")
    plt.legend()
    plt.grid()
    plt.savefig(pic)
    plt.show()






