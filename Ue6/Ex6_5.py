# import math as m
# import time
# import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

def plotError_loglog(fig, ax, X, Y, exact=0, filepath="/Users/peterholzner/Code/NumComp/newplotfile", style="x-", label = "no label", plot_title='loglog plot', horder=0, save=0, debug=0):
    # n = len(X)
    ax.loglog(X, list(map(lambda x: abs(x - exact), Y)), style, label=label)
    order = int(horder)
    while(order > 0):
        ax.loglog(X,list(map(lambda h: h**order, X)), '-', label="h^"+str(order))
        order -= 1
    ax.legend()
    ax.set(xlabel='Knots h_i', ylabel='error', title=plot_title)
    ax.autoscale
    # i = ax.get_ylim()
    # if i[0] < 10**(-17):
    #     ax.set_ylim(bottom=10**(-17))
    #     print(i)
    # ax.set(ylim=(10**-16, 1))
    ax.grid(b= True)
    if save>0: fig.savefig(filepath)
    if debug>0: plt.show()

def plotError_semilogy(fig, ax, X, Y, exact=0, filepath="/Users/peterholzner/Code/NumComp/newplotfile", style="x-",label = "no label", plot_title='loglog plot', horder=0, save=0, debug=0):
    # n = len(X)
    ax.semilogy(X, list(map(lambda x: abs(x - exact), Y)), style, label=label)
    order = int(horder)
    while(order > 0):
        ax.semilogy(X,list(map(lambda h: h**order, X)), '-', label="h^"+str(order))
        order -= 1
    ax.legend()
    ax.set(xlabel='Knots h_i', ylabel='error', title=plot_title)
    ax.autoscale
    i = ax.get_ylim()
    if i[0] < 10**(-17):
        ax.set_ylim(bottom=10**(-17))
    #     print(i)
    # ax.set(ylim=(10**-16, 1))
    ax.grid(b= True)
    if save>0: fig.savefig(filepath)
    if debug>0: plt.show()

def aitken(a,b,c):
    return a - ((b-a)**2) / (a-2*b+c)

if __name__ == "__main__":
    pngpath = "/Users/peterholzner/Code/NumComp/Ue6/Ex6_5_error.pdf"
    kMax = 30
    iMax = int(kMax/3)
    eps = 1e-8
    safety = 2
    iConvCheck = 4
    print("kMax = ", kMax)
    Ks = [k for k in range(0,kMax+1)]
    U = np.zeros_like(Ks, dtype=float)
    A = np.zeros_like(U)
    U[1] = 2
    A[0] = 2
    A[1] = 2
    iFirst = True
    kFirst = True
    for k in Ks[1:-1]:
        U[k+1] = (2**k) * np.sqrt(2*(1 - np.sqrt(1 - ((0.5**k) * U[k])**2)))
        if abs(U[k+1] - U[k]) < eps: # here with criterion that difference to previous element is < eps=1e-8
            k_eps = k+1
        if (U[k+1]-U[k]) * (U[k] - U[k-1]) < 0 and kFirst: # here with criterion is sign change of difference to previus element
            k_sign = k
            kFirst = False
    i = 0
    A[i] = aitken(U[i+1], U[i+2], U[i+3])
    for i in Ks[1:-2]:
        A[i] = aitken(U[i], U[i+1], U[i+2])
        if abs(A[i] - A[i-1]) < eps: # here with criterion that difference to previous element is < eps=1e-8
            i_eps = i
        if i > 1: # here with criterion is sign change of difference to previus element
            if ((A[i]-A[i-1]) * (A[i-1] - A[i-2]) < 0 or abs(A[i]-A[i-1])>abs(A[i-1] - A[i-2])*safety) and iFirst and i>iConvCheck:
                i_sign = i - 1  
                iFirst = False

    exact = np.pi
    errorsU = np.absolute(U - exact)
    errorsA = np.absolute(A - exact)
    print("U[kMax=", kMax, "] = ", U[kMax])
    print("Exact result = pi = ", exact)
    print("Error of last member of sequence: ", errorsU[kMax])
    print("Error of last member of aitken sequence: ", errorsA[iMax])
    print("-"*10)
    print("with eps = ", eps)
    print("Original sequence u[k] stops at u[k=", k_eps,"]= ", U[k_eps])
    print("Aitken sequence a[i] stops at a[i=", i_eps,"]= ", A[i_eps])
    print("-"*10)
    print("with sign change criterion:")
    print("Original sequence u[k] stops at u[k=", k_sign,"]= ", U[k_sign])
    print("Aitken sequence a[i] stops at a[i=", i_sign,"]= ", A[i_sign])
    fig, ax = plt.subplots()
    plotError_semilogy(fig, ax, Ks, errorsU, label="u[k]",filepath=pngpath, save=1)
    plotError_semilogy(fig, ax, Ks, errorsA, label="aitken seq a[i]",filepath=pngpath, save=1)
    plt.close()