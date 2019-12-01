import math as m
import time
import scipy.integrate as integrate
import matplotlib.pyplot as plt
#import numpy as np

const1 = 1/3

def f1(x):
    return x**2

def F1(x):
    return (x**3)/3

def f2(x):
    # second derivatie doesnt exist --> trapezoidal is exact
    return abs(x)
    

def F2(x):
    if x>=0:
        return 0.5*f1(x)
    else:
        return -0.5*f1(x)

def f3(x):
    # second deriv is smth like exp(x) again
    if x< 1/3:
        return 0.5*m.exp(x)
    else:
        return m.exp(x)

def F3(x):
    return f3(x)

def integrateTrap(N, f=f1, a=-1, b=1):
    h = (b-a)/N
    t_list = list(map(lambda x: a + x*h, list(range(1,int(N)))))
    return h*(0.5*(f(a)+f(b)) + sum(map(f, t_list)))

def plotError_loglog(X, Y, X2=[], Y2=[], exact=0, filepath="/Users/peterholzner/Code/NumComp/newplotfile", labels = ['plot1', 'plot2'], plot_title='loglog plot'):
    n = len(X2)
    fig, ax = plt.subplots()
    
    ax.loglog(X, list(map(lambda x: abs(x - exact), Y)), 'x-', label=labels[0])
    if n:
        ax.loglog(X2, list(map(lambda x: abs(x-  exact), Y2)), 's-', label=labels[1])
    ax.legend()
    ax.set(xlabel='Knots h_i', ylabel='error', title=plot_title)
    # i = ax.get_ylim()
    # if i[0] < 10**(-17):
    #     ax.set_ylim(bottom=10**(-17))
    #     print(i)
    # ax.set(ylim=(10**-6, 1))
    ax.grid()
    # Save fig into file (BEWARE: absolute path!)
    fig.savefig(filepath)
    # plt.show()

if __name__ == "__main__":
    a= -1
    b = 1
    i_max = 20
    iters = list(range(0,i_max+1))
    h_list = list(map(lambda x: 2**-x, iters))
    N_list = list(map(lambda x: (b-a)/x, h_list))
    for f,F in zip((f1, f2, f3),(F1, F2, F3)):
        # Integration by hand
        
        if f==f3:
            eps = 0.00001
            #exactVal = F(b) - F(const1+eps) + F(const1-eps) - F(a)
            exactVal = 0.5*m.exp(1/3)-0.5*m.exp(-1)+m.exp(1)-m.exp(1/3)
        else:
            exactVal = F(b) - F(a)

        # Scipy Integration
        start = time.time()
        result_scipy = integrate.quad(f, a, b)
        rt_scipy = time.time() - start
        
        # Custom routine
        start = time.time()
        result_trap = integrateTrap(N_list[i_max], f, a, b)
        rt_trap = time.time() - start

        results = []
        for n in N_list:
            results.append(integrateTrap(n, f=f))
        print(results)
        # Output
        if f==f1:
            print("f(x)=x^2 --> F(x)=(x^3)/3 +C")
            name = "f1(x)=x^2"
            path = "/Users/peterholzner/Code/NumComp/Ue3/f1"
        elif f==f2:
            print("f(x)=abs(x) --> F(x)=signum(x)*(x^2)/2 +C")
            name = "f2(x)=abs(x)"
            path = "/Users/peterholzner/Code/NumComp/Ue3/f2"
            
        elif f==f3:
            # eps = 10^-8
            # exactVal = F(b) - F(const1+eps) + F(const1-eps) - F(a)
            print("f(x)=c*e^x (c=1 if x>=1/3 or c=0.5 else) --> F(x)=f(x)")
            name = "f3(x)=c*e^x (c=1 if x>=1/3 or c=0.5 else)"
            path = "/Users/peterholzner/Code/NumComp/Ue3/f3"
            # print(results)
        print("Exact value: ", exactVal)
        print("Result from the scipy integration routine quad:")
        print(result_scipy)
        print("Runtime: ", rt_scipy)
        print("Result from my custom trapezoid integration routine with (N=", N_list[i_max], "): ")
        print(result_trap)
        print("Runtime: ", rt_trap)
        #for i in list(map(lambda x: x*2, range(0,11))): print(results[i], "\n")
        plotError_loglog(h_list, results, exact= exactVal, filepath= path,plot_title=name)
        if f==f1 or f==f2:
            print("______________________________")