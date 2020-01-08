import math as m
# import time
# import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
import os

def f0(x):
    return x**0

def f1(x):
    return x

def f2(x):
    return x**2

def fb(x):
    return np.log(x)*np.power(x,0.1)

def f_fitlaw(x, C, b):
    return C*np.exp(-b*x)

def f_54a_ln(x): # log(x) is bad for small x...so my substitution isn't suited for numerical computation (it's not stable?)
    return np.log(1-np.log(x))/(x*np.power(1-np.log(x),np.pi))

def f_54a_1_u(x):
    return (-1)*np.log(x) * np.power(x, np.pi-2)

def f_54b(x):
    return x/np.exp(x*(np.pi - 1))

def comp_gauss(func, n, L, q, a=0, b=1, info=0):
    """n...points per subintervall\n
    L...number of subintervalls
    q...quadrature base
    """
    if q>=1: 
        print("Please select a q<1!\nReturning 0!")
        return 0
    Is = b*np.array([0]+[np.power(q,i) for i in np.arange(L-1, -1, -1)])
    lenIs = np.array([p2-p1 for p1, p2 in zip(Is[:-1], Is[1:])])
    Xg = np.zeros((L,n))
    Wg = np.zeros((1,n))
    xp, Wg = np.polynomial.legendre.leggauss(n)
    # xp = np.array(xp)
    for x, a, b, len in zip(Xg,Is[:-1], Is[1:], lenIs):
        x += 0.5*(len*xp+a+b)  
    if info:
        print("Intervalls:\n", Is)
        print("Length of intervalls:\n", lenIs)
        print("xg:\n", Xg)
        print("wg:\n", Wg)
    res = 0.5 * np.dot(lenIs, np.dot(func(Xg), Wg))
    return res

def cedrik_comp_gauss(f, n, L, q, a=0, b=1):
    ret = 0
    b = b * (q ** (L - 1))
    x, w = np.polynomial.legendre.leggauss(n)
    for i in range(L):
        temp = i
        temp += 1
        ret += (b - a) / 2 * sum(w * f((b - a)/2 * x + (a + b)/2))
        a = b
        b = b / q
    return ret

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
    ax.grid()
    if save>0: fig.savefig(filepath)
    if debug>0: plt.show()

def plotError_semilogy(fig, ax, X, Y, exact=0, filepath="/Users/peterholzner/Code/NumComp/newplotfile", style="x-",label = "no label", plot_title='loglog plot', horder=0, save=0, debug=0):
    # n = len(X)
    ax.semilogy(X, np.absolute(Y - exact), style, label=label)
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
    ax.grid()
    if save>0: fig.savefig(filepath)
    if debug>0: plt.show()

if __name__ == "__main__":
    scriptpath = os.path.dirname(__file__)
    filepath = os.path.join(scriptpath, "Ex5_1_error")
    fig, ax = plt.subplots()
    fitfunc = np.polynomial.polynomial.Polynomial
    #####################
    n = 3
    L = 5
    q = 0.5
    print("Testing my composite gauss function for 3 test functions:\nf0(x)=1,\tf1(x)=x,\tf2(x)=x^2")
    exact = [1, 0.5, 1/3]
    res = [comp_gauss(func, n, L, q, a=0, b=1) for func in (f0, f1, f2)]
    # resDL = [comp_gauss(func, n, L, q, a=0, b=1)/L for func in (f0, f1, f2)]
    print("Exact:\n", exact)
    print("Result\n", res)
    # print("Result/L\n", resDL)
    #---------------5.1.b+c-----------------
    print("New function!\nf(x) = log(x)/x^π")
    exact = -1/(1.1**2)
    n_max = 20

    Ns = np.arange(1, n_max+1)
    Results = np.zeros((3,n_max))
    Errors = np.zeros_like(Results)
    Qs = [0.5, 0.15, 0.05]
    nr_qs = len(Qs)
    for q, i in zip(Qs, np.arange(0,nr_qs)):
        temp = np.zeros(n_max)
        temp = np.array([comp_gauss(fb, n, n, q, a=0, b=1) for n in Ns])
        error_temp = np.absolute(temp - exact)
        plotError_semilogy(fig, ax, Ns, error_temp, exact=0, filepath=filepath,label="q="+str(q), plot_title="comp_gauss convergence", save=0)
        Results[i,:] = temp
        Errors[i,:] = error_temp

    for i, errs in zip(np.arange(0,nr_qs), Errors):
        Errors_log = np.log(errs[errs!=0])
        if i == 1:
            Ns1 = Ns[Ns != 18]
            coeffs = fitfunc.fit(Ns1, Errors_log, 1)
            c_poly = np.polyfit(Ns1, Errors_log,1)
            print("Coeffs: ", c_poly)
        else:    
            coeffs = fitfunc.fit(Ns, Errors_log, 1)
            c_poly = np.polyfit(Ns, Errors_log, 1)
            print("Coeffs: ", c_poly)
        print(coeffs)
        # C_coef = np.exp(coeffs.coef[1])
        # b_coef = -coeffs.coef[0]
        C_coef = np.exp(c_poly[1])
        b_coef = -c_poly[0]
        error_fi = f_fitlaw(Ns, C=C_coef, b=b_coef)
        plotError_semilogy(fig, ax, Ns, error_fi, exact=0, filepath=filepath, style='--',label=str(i), plot_title="comp_gauss convergence", save=1)
        
    # q=0.15 seems to be the best choice of the 3
    # Errors = np.zeros_like(Results)
    # Errors = Results - exact
    # Errors = np.absolute(Errors)
        # "C="+str(C_coef)+"\nb="+str(b_coef)
    # fit to C*e^-b*n
    # log(...) = log(C) + (-b)*log(e)*n = log(C) + (-b)*n //log=ln
    # ==> a0 = log(C), a1 = -b
    # print(Errors)
    
    print("\n---------Ex5.4---------\n")
    exact = -1/(1.1**2)
    n = 5
    L = 5
    q = 0.15
    print("Testing my composite gauss function for 3 test functions:\nf0(x)=f_5_4_a")
    exact = 1/(np.power(np.pi, 2) - 2*np.pi + 1)
    resa1 = comp_gauss(f_54a_ln, n, L, q, a=0, b=1)
    resa2 = comp_gauss(f_54a_1_u, n, L, q, a=0, b=1)
    n = 7
    L = 10
    q = 0.1
    res2 = comp_gauss(f_54b, n, L, q, a=0, b=L)
    print("Exact:\n", exact)
    print("Result for f_54a_ln\n", resa1)
    print("Result for f_54a_1-u\n", resa2)
    print("Result2\n", res2)

    q = 0.15
    fig, ax = plt.subplots()
    filepath = os.path.join(scriptpath, "Ex5_4_error")
    for func, end, lab in zip((f_54a_ln, f_54a_1_u, f_54b), (1, 1, 100) ,["u=e^(1-x)", "u=1/x", "F(y)"]):
        temp = np.zeros(n_max)
        temp = np.array([comp_gauss(func, n, n, q, a=0, b=end) for n in Ns])
        error_temp = np.absolute(temp - exact)
        plotError_semilogy(fig, ax, Ns, error_temp, exact=0, filepath=filepath,label=lab, plot_title="comp_gauss convergence", save=1)

    

    
    
            