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

if __name__ == "__main__":
    pngpath = "/Users/peterholzner/Code/NumComp/Ue6/Ex6_4_error"
    print("Please input a Kmax > 0!")
    kMax = int(input())
    if kMax < 0:
        print("kMax = ", kMax, " will not be used! Only positive kMax allowed!")
        kMax = 20
        print("Using kMax = ", kMax, " instead.")
    Ks = [k for k in range(0,kMax+1)]
    U = np.zeros_like(Ks, dtype=float)
    U[1] = 2
    for k in Ks[1:-1]:
        U[k+1] = (2**k) * np.sqrt(2*(1 - np.sqrt(1 - ((0.5**k) * U[k])**2)))

    exact = np.pi
    errors = np.absolute(U - exact)
    print("U[kMax=", kMax, "] = ", U[kMax])
    print("Exact result = pi = ", exact)
    print("Error of last member of sequence: ", errors[kMax])
    print(U)
    fig, ax = plt.subplots()
    plotError_semilogy(fig, ax, Ks, errors, filepath=pngpath, save=1)
    plt.show()
    plt.close()