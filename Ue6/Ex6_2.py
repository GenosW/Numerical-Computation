# import math as m
# import time
# import scipy.integrate as integrate
# import numpy as np
import matplotlib.pyplot as plt

def f1(x,y):
    return x**2

def f2(x, y):
    if x < y: return 0
    else: return 1

def lenInt(I):
    return I[1] - I[0]

def divideInt(I):
    mid = (I[0]+I[1])*0.5
    I1 = [I[0], mid]
    I2 = [mid, I[1]]
    return I1, I2

def adaptRect(f, Ix, Iy, tau, hmin):
    # Check if minimum partioning has been reached
    Msq = Midpoint(f, Ix, Iy)
    if lenInt(Ix) < hmin or lenInt(Iy) < hmin: return Msq
    # Finer partitioning requested
    Ix1, Ix2 = divideInt(Ix)
    Iy1, Iy2 = divideInt(Iy)
    Mm = sum([Midpoint(f, i, j) for i in (Ix1, Ix2) for j in (Iy1, Iy2)])
    # Check if finer partioning improved the result
    if abs(Msq - Mm) <= tau: return Msq # Already good enough
    else: # Use finer partioning
        newtau = tau*0.25
        res = sum([adaptRect(f, i, j, newtau, hmin) for i in (Ix1, Ix2) for j in (Iy1, Iy2)])
        return res

def Midpoint(f, Ix, Iy):
    return (Ix[1]-Ix[0])*(Iy[1]-Iy[0])*f(0.5*(Ix[0]+Ix[1]), 0.5*(Iy[0]+Iy[1]))

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
    a = 0
    b = 1
    c = 0
    d = 1
    Ix = [a,b]
    Iy = [c,d]
    i_max = 15
    taus = list(map(lambda i: 0.5**i, range(0,i_max+1)))
    pngpath = "/Users/peterholzner/Code/NumComp/Ue6/Ex6_2_error"
    
    res1 = [adaptRect(f1, Ix, Iy, tau, tau) for tau in taus] 
    exact1 = 1/3
    res2 = [adaptRect(f2, Ix, Iy, tau, tau) for tau in taus]
    exact2 = 0.5

    fig, ax = plt.subplots()
    print("Integral over the unit square in 2D of f1(x,y)=x^2")
    print("Result of adaptRect: ", res1)
    print("Exact value: ", exact1)
    plotError_loglog(fig, ax, taus, res1, exact=exact1, filepath=pngpath, label="f1(x)", save=0)

    print("Integral over the unit square in 2D of f2(x,y)= (0 if x<y) (1 if x>y)")
    print("Result of adaptRect: ", res2)
    print("Exact value: ", exact2)
    plotError_loglog(fig, ax, taus, res2, exact=exact2, filepath=pngpath, label="f2(x)", plot_title="Convergence of adaptive Integration over the unit square",horder=1, save=1)
    plt.close()