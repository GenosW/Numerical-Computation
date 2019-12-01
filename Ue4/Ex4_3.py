import math as m
import time
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np

def f1(x):
    return x**0.2

def f2(x):
    return x**10

def f3(x):
    return x**2

def F(x, a=0.2):
    return x**(1+a)/(1+a)

def integrateTrap(h, f=m.exp, a=0, b=1):
    # h = (b-a)/N
    N = (b-a)/h
    t_list = list(map(lambda x: a + x*h, list(range(1,int(N)))))
    return h*(0.5*(f(a)+f(b)) + sum(map(f, t_list)))

def plotError_loglog(fig, ax, X, Y, exact=0, filepath="/Users/peterholzner/Code/NumComp/newplotfile", label = "no label", plot_title='loglog plot', debug=0):
    # n = len(X)
    ax.loglog(X, list(map(lambda x: abs(x - exact), Y)), 'x-', label=label)
    ax.legend()
    ax.set(xlabel='Knots h_i', ylabel='error', title=plot_title)
    # i = ax.get_ylim()
    # if i[0] < 10**(-17):
    #     ax.set_ylim(bottom=10**(-17))
    #     print(i)
    ax.set(ylim=(10**-16, 1))
    ax.grid()
    # Save fig into file (BEWARE: absolute path!)
    fig.savefig(filepath)
    if debug:
        plt.show()

def neville(x, f, point=0, n=-1, retVal = "scheme"):
    """ Neville scheme:\n
    x...        vector of knots of length n+1,\n
    f...        vector of data values of length n+1,\n
    point...    the point(knot) to evaluate the neville scheme at,\n
    n...        size of the neville scheme, calculated automatically if left at default,\n
    retVal...   determines what is given as return value (scheme or scheme + value(point)).
    """
    if n<0:
        x = np.array(x)
        n = x.size
    q = np.zeros((n, n))
    q[:,0] = np.copy(f)

    for m in range(1, n):
        for i in range(0, n-m):
            q[i,m] = (point - x[i])*q[i+1,m-1] - (point - x[i+m])*q[i,m-1]
            q[i,m] /= (x[i+m] - x[i])
    if retVal=="point":
        result = q[i,m]
    elif retVal=="scheme":
        result = q
    else:
        result = q
    return result

if __name__ == "__main__":
    a= 0
    b = 1
    i_max = 10
    alpha = [0.2, 10, 2]
    iters = list(range(0,i_max+1))
    h_list = list(map(lambda x: 0.5**x, iters))
    h2 = list(map(lambda x: 0.5**(2*x), iters))
    pngpaths = ["/Users/peterholzner/Code/NumComp/Ue4/Ex4_3_error_f" + str(i) for i in range(1,4)]

    for func, alpha, path in zip((f1, f2, f3), (0.2, 10, 2), pngpaths):
        fi = [integrateTrap(h, f=func) for h in h_list]
        exact = F(b, a=alpha) - F(a, a=alpha)
        res = neville(h2, fi)

        fig, ax = plt.subplots()
        i_liste = [0,1,2]
        print(exact)
        for i in i_liste:
            biggestError = abs(exact - res[0,i])
            smallestError = abs(exact - res[i_max-i,i])
            if biggestError < 10**(-16):
                print("For x**",alpha,"\n the neville column m=", i, "is exact!")
            if func == f3 and i == 2:
                print(biggestError)
            if func == f3 and i==1:
                print(res[0:i_max-i,i])
            plotError_loglog(fig, ax, h_list[0:i_max-i], res[0:i_max-i,i], exact=exact, filepath=path, label="m="+str(i))
        plt.close()