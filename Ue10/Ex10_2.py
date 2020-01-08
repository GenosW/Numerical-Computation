import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
# import random

def f(x):
    assert(isinstance(x,np.ndarray))
    return np.array([3*x[0] - np.cos(x[1]*x[2]) - 1.5,\
                    4*x[0]*x[0] - 625*x[1]*x[1] + 2*x[2] - 1,\
                    20*x[2] + np.exp(-x[0]*x[1]) + 9], dtype=float)


def df(x):
    assert(isinstance(x,np.ndarray))
    return np.array( [[3, x[2]*np.sin(x[1]*x[2]), x[1]*np.sin(x[1]*x[2])],\
                    [8*x[0], -1250*x[1], 2],\
                    [-x[1]*np.exp(-x[0]*x[1]), -x[0]*np.exp(-x[0]*x[1]), 20]], dtype = float)

def netwonXD_inv(xn,f,df):
    if np.linalg.det(df(xn)):
        inv = np.linalg.inv(df(xn))
        fxn = f(xn)
        # print("inv:", inv)
        # print("fxn", fxn)
        return xn - np.dot(inv,fxn)
    else:
        print("df is noninvertible...returning input point x_n!")
        return xn

def netwonXD(xn,f,df):
    if np.linalg.det(df(xn)):
        res = f(xn)
        corr = np.linalg.solve(df(xn),res)
        return xn - corr, res
    else:
        print("df is noninvertible->cant be solved...returning input point x_n!")
        return xn

def vecNorm(x):
    return np.sqrt(np.sum(x*x))

if __name__ == "__main__":
    x = np.array([1,1,1], dtype=float)
    print("Starting guess x0: ",x)
    print("norm(x0): ",np.linalg.norm(f(x)))
    print("-"*30)
    print("Iterating...please hold the line!")
    eps = 1e-16
    error = [np.linalg.norm(f(x))]
    iterations = 100
    for i in list(range(0,100)):
        xn1, res = netwonXD(x,f,df)
        error.append(np.linalg.norm(res))
        if np.allclose(xn1,x,rtol=eps):
            # print("x:\n",x)
            # print("xn1:\n",xn1)
            print("Break condition: <xn, xnp1 close>")
            iterations = i+1
            np.copyto(x, xn1)
            break
        np.copyto(x, xn1)
    print("Arrived at solution after ",iterations," iterations.")
    print("final:\n",x)
    print("final error:\n",np.linalg.norm(f(x)))
    print("-"*30)
    print("Trying different break condition...please hold the line")
    x = np.array([1,1,1], dtype=float)
    error = [np.linalg.norm(f(x))]
    iterations = 100
    for i in list(range(0,100)):
        xn1, res = netwonXD(x,f,df)
        error.append(np.linalg.norm(res))
        if np.abs(np.linalg.norm(x) - np.linalg.norm(xn1)) < eps:
            # print("x:\n",x)
            # print("xn1:\n",xn1)
            print("Break condition: <norms converged>")
            iterations = i+1
            break
        np.copyto(x, xn1)
    print("Arrived at solution after ",iterations," iterations.")
    print("final:\n",x)
    print("final error:\n",np.linalg.norm(f(x)))
    x2 = np.power(1e-2,np.arange(0,iterations+1))
    plt.semilogy(np.arange(0,iterations+1),error,"rx-")
    #plt.semilogy(np.arange(0,iterations+1),x2,"b-")
    plt.grid()
    scriptpath = os.path.dirname(__file__)
    plt.savefig(os.path.join(scriptpath,"Ex10_2.png")) 
    
        
        
    