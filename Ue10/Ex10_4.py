import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
# import random

def netwonXD(A,xn,f,df,n,h):
    if np.linalg.det(df(xn,n,h)):
        res = f(xn,n,h,A)
        corr = np.linalg.solve(df(xn,n,h),res)
        xNew = xn - corr
        return xNew#, f(xNew,n,h,A)
    else:
        print("df is noninvertible->cant be solved...returning input point x_n!")
        return xn, 0*xn

def f(u,n,h,A):
    ret = np.dot(A,u) + np.power(u,3) - np.ones(n,dtype=float)
    return ret

def df(u,n,h):
    ret = np.diagflat(3*np.power(u,2) +2/h**2) - 1/h**2 * (np.eye(n,n,k=1) + np.eye(n,n,k=-1))
    return ret

if __name__ == "__main__":
    # n = 10
    n = int(input("n = "))
    print("n = ",n," selected!")
    h = 1/(n-1)
    A = 1/h**2 * (2*np.eye(n,n,dtype=float) - np.eye(n,n,k=1,dtype=float) - np.eye(n,n,k=-1,dtype=float))
    b = np.zeros(n,dtype=float)
    ########
    print("-"*30)
    print("Now Newton: ")
    #x = np.linalg.solve(A,b)
    x = np.zeros(n,dtype=float)
    # print("Starting guess x0: ",x)
    print("norm(g(x0)): ",np.linalg.norm(f(x,n,h,A)))
    print("-"*30)
    print("Iterating...please hold the line!")
    eps = 1e-15
    error = []
    iterations = 100
    for i in list(range(0,iterations)):
        xn1 = netwonXD(A,x,f,df,n,h)
        error.append(np.linalg.norm(xn1 - x))
        if error[-1] < eps:# or np.allclose(x,xn1,rtol=eps):
            # print("Break condition: <xn, xnp1 close>")
            iterations = i+1
            np.copyto(x,xn1)
            break
        np.copyto(x, xn1)
    close = np.allclose(f(x,n,h,A),b,rtol=eps)
    print("Is it closed? ",close)
    print("Iterations: ", iterations)
    # print("final form of x: ",x)
    # print(error)
    plt.semilogy(np.arange(0,iterations,step=1),error,"bx-",label="Newton")
    plt.grid()
    plt.legend()
    # plt.show()
    scriptpath = os.path.dirname(__file__)
    plt.savefig(os.path.join(scriptpath,"Ex10_4.png")) 

