import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
# import random

class simple:
    def __init__(self, A, b, eps, f):
        self.A = A
        self.b = b
        self.eps = eps
        self.f = f
        self.x0 = self.getStart()
        self.x = self.x0
        self.iterations = 0
        self.xList = [self.x]
        self.errors = []
        self.estErrors = []

    def getStart(self):
        return np.linalg.solve(self.A,self.b)

    def rhs(self,x):
        return self.b + self.eps*self.f(x)

    def estimateError(self,i):
        return np.linalg.norm(self.xList[-i]-self.xList[-i-1])

    def iterate(self,maxIterations,rtol):
        for i in list(range(0,maxIterations)):
            rhsBef = self.rhs(self.x)
            self.x = np.linalg.solve(self.A, rhsBef)
            self.xList.append(self.x)
            est = self.estimateError(1)
            # if np.allclose(np.dot(self.A,self.x),self.rhs(self.x),rtol=rtol):
            if (est < rtol): # gives error of |xn+1 - xn|_2
                self.iterations = i
                break
            self.estErrors.append(est)
        return self.x

    def close(self,rtol):
        return np.allclose(np.dot(self.A,self.x),self.rhs(self.x),rtol=rtol)

    def calcErrors(self):
        for xi in self.xList[:-2]:
            self.errors.append(np.linalg.norm(self.x - xi))

def f(x):
    return np.array([(x[0]-x[1])**2, 0], dtype=float)

def g(x):
    return np.array([2*x[0]+x[1]-0.01*(x[0]-x[1])**2,\
                    x[0]+2*x[1]-3], dtype=float)

def dg(x):
    diff = x[0]-x[1]
    return np.array( [[2*(1-0.01*diff), 1+0.02*diff],\
                        [1, 2]], dtype = float)

def netwonXD(xn,f,df):
    if np.linalg.det(df(xn)):
        res = f(xn)
        corr = np.linalg.solve(df(xn),res)
        xNew = xn - corr
        return xNew, f(xNew)
    else:
        print("df is noninvertible->cant be solved...returning input point x_n!")
        return xn, 0*xn

if __name__ == "__main__":
    print("Simple first:")
    A = np.array([[2,1],[1,2]], dtype=float)
    b = np.array([0,3], dtype=float)
    eps = 1e-2
    maxIterations = 100
    eq = simple(A,b,eps,f)
    xFin = eq.iterate(maxIterations,1e-16)
    print("final:\n",xFin)
    print(eq.close(1e-16))
    print(np.linalg.norm(np.dot(A,xFin) - eq.rhs(xFin)))
    eq.calcErrors()
    print("Errors: ")
    print(eq.estErrors)
    ##########################
    print("-"*30)
    print("Now Newton: ")
    x = np.linalg.solve(A,b)
    print("Starting guess x0: ",x)
    print("norm(g(x0)): ",np.linalg.norm(g(x)))
    print("-"*30)
    print("Iterating...please hold the line!")
    eps = 1e-16
    error = [np.linalg.norm(g(x))]
    iterations = 100
    for i in list(range(0,100)):
        xn1, res = netwonXD(x,g,dg)
        error.append(np.linalg.norm(res))
        if np.linalg.norm(res) < eps or np.allclose(x,xn1,rtol=eps):
            print("Break condition: <xn, xnp1 close>")
            iterations = i+1
            np.copyto(x,xn1)
            break
        np.copyto(x, xn1)
    print("Arrived at solution after ",iterations," iterations.")
    print("final:\n",x)
    print("final error:\n",np.linalg.norm(np.dot(A,x) - g(x)))
    #x2 = np.power(1e-2,np.arange(0,iterations+1))
    plt.semilogy(np.arange(0,iterations+1),error,"bx-",label="Newton")
    plt.semilogy(np.arange(0,len(eq.errors)),eq.errors,"go-",label="Simple_actual error")
    plt.semilogy(np.arange(0,len(eq.estErrors)),eq.estErrors,"rx--",label="Simple_estimated error")
    # plt.semilogy(np.arange(0,iterations+1),x2,"b-", label="order x^2 ???")
    plt.grid()
    plt.legend()
    # plt.show()
    scriptpath = os.path.dirname(__file__)
    plt.savefig(os.path.join(scriptpath,"Ex10_3.png")) 
