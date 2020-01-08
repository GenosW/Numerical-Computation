import numpy as np
import matplotlib.pyplot as plt
import os
# import random

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

def newton(x, f, df):
    return x - f(x)/df(x)

def iterateNewton(x0, f, df, maxIteration, eps, bound=1000, debug=0):
    x = newton(x0,f,df)
    xList = [x0,x]
    iterations = 0
    for i in list(range(0,maxIteration)):
        if debug: print("Iteration: ", i+1,"\nx= ", x)
        x = newton(x,f,df)
        xList.append(x)
        if (abs(f(x)) < eps):
            iterations = i+1
            break
        if (abs(x-x0) >= bound):
            print("Out of bound! Terminating Newton-Raphson method.")
            iterations = i+1
            break
    return xList, iterations


def f1(x):
    return np.power(x,2)

def df1(x):
    return 2*x

def f2(x):
    return np.exp(x) - 2

def df2(x):
    return np.exp(x)

def f3(x):
    return np.power(np.abs(x), 1.5)

def df3(x):
    return 1.5*np.sqrt(np.abs(x))*np.sign(x)

def f4(x):
    return 1/x - 1

def df4(x):
    return -1/np.power(x,2)

if __name__ == "__main__":
    eps = 1e-16
    x0 = 0.5
    x02 = 2.1
    # sigPrev = np.sign(f(x))
    print("f1")
    xExact1 = 0
    xList1, iterations1 = iterateNewton(x0, f1, df1, 50, eps)
    print("x_root = ", xList1[-1])
    print("found after ",iterations1," iterations")
    print("with an accuracy of", abs(xExact1 - xList1[-1]))
    print("f2")
    xExact2 = np.log(2)
    xList2, iterations2 = iterateNewton(x0, f2, df2, 50, eps)
    print("x_root = ", xList2[-1])
    print("found after ",iterations2," iterations")
    print("with an accuracy of", abs(xExact2 - xList2[-1]))
    print("f3")
    xExact3 = 0
    xList3, iterations3 = iterateNewton(x0, f3, df3, 50, eps)
    print("x_root = ", xList3[-1])
    print("found after ",iterations3," iterations")
    print("with an accuracy of", abs(xExact3 - xList3[-1]))
    print("f4")
    xExact4 = 1
    xList4, iterations4 = iterateNewton(x02, f4, df4, 50, eps, bound=10)
    print("x_root = ", xList4[-1])
    print("found after ",iterations4," iterations")
    print("with an accuracy of", abs(xExact4 - xList4[-1]))
    print(xList4)

    fig, ax = plt.subplots()
    # plotError_semilogy(fig, ax, list(range(0,iterations1+2)), xList1, xExact1, save=1)
    ax.semilogy(list(range(0,iterations1+2)), np.abs(np.array(xList1) - xExact1), '-x', label="x²")
    ax.semilogy(list(range(0,iterations2+2)), np.abs(np.array(xList2) - xExact2), '-x', label="e^x -1")
    ax.semilogy(list(range(0,iterations3+2)), np.abs(np.array(xList3) - xExact3), '-x', label="|x|^3/2")
    ax.semilogy(list(range(0,iterations4+2)),np.abs( np.array(xList4) - xExact4), '-x', label="x^-1 - 1")
    # ax.semilogy(list(range(0,max(iterations1,iterations2,iterations3)+2)), )
    ax.grid()
    ax.legend()
    ax.set(title="Yo", xlabel="iterations", ylabel="error of x_root")
    # plt.show()
    
    # print(os.path.join(scriptpath,"Ex9_2.png")
    scriptpath = os.path.dirname(__file__)
    fig.savefig(os.path.join(scriptpath,"Ex9_2.png"))    