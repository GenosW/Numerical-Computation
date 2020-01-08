import math as m
import time
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np

alpha = 1/3

def f(x):
    if x < alpha: return 0.5*m.exp(x)
    else: return m.exp(x)

def adapt(f, a, b, tau, hmin):
    Sab = Simpson(f, a, b)
    if (b-a) < hmin: return Sab
    m = (a+b)*0.5
    Samb = Simpson(f, m, b) + Simpson(f, a, m)
    if abs(Sab - Samb) <= tau: return Sab
    else: 
        newtau = tau*0.5
        res = adapt(f, a, m, newtau, hmin) + adapt(f, m, b, newtau, hmin)
    return res

def Simpson(f, a, b):
    return (b-a)/6 * (f(a) + 4*f((a+b)*0.5) + f(b))

def plotError_loglog(fig, ax, X, Y, exact=0, filepath="/Users/peterholzner/Code/NumComp/newplotfile", label = "no label", plot_title='loglog plot', horder=0, debug=0):
    # n = len(X)
    ax.loglog(X, list(map(lambda x: abs(x - exact), Y)), 'x-', label=label)
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
    ax.grid()
    # Save fig into file (BEWARE: absolute path!)
    fig.savefig(filepath)
    if debug:
        plt.show()

# def neville(x, f, point=0, n=-1, retVal = "scheme"):
#     """ Neville scheme:\n
#     x...        vector of knots of length n+1,\n
#     f...        vector of data values of length n+1,\n
#     point...    the point(knot) to evaluate the neville scheme at,\n
#     n...        size of the neville scheme, calculated automatically if left at default,\n
#     retVal...   determines what is given as return value (scheme or scheme + value(point)).
#     """
#     if n<0:
#         x = np.array(x)
#         n = x.size
#     q = np.zeros((n, n))
#     q[:,0] = np.copy(f)

#     for m in range(1, n):
#         for i in range(0, n-m):
#             q[i,m] = (point - x[i])*q[i+1,m-1] - (point - x[i+m])*q[i,m-1]
#             q[i,m] /= (x[i+m] - x[i])
#     if retVal=="point":
#         result = q[i,m]
#     elif retVal=="scheme":
#         result = q
#     else:
#         result = q
#     return result

if __name__ == "__main__":
    a= 0
    b = 1
    i_max = 30
    iters = list(range(0,i_max+1))
    hmins = list(map(lambda j: 0.5**j, iters))
    pngpath = "/Users/peterholzner/Code/NumComp/Ue4/Ex4_4_error"
    
    I = [adapt(f, a, b, h, h) for h in hmins] 
    exact = m.exp(1) - 0.5*(1+m.exp(1/3))
    print(I)

    fig, ax = plt.subplots()
    print("Exact value: ",exact)
    plotError_loglog(fig, ax, hmins, I, exact=exact, filepath=pngpath, label="f(x)", plot_title="Convergence of adaptive Simpson",horder=1)
    plt.close()