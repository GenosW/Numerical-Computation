import os
import numpy as np
import matplotlib.pyplot as plt 

def F(x):
    return 2 - np.power(x,2) - np.exp(x)

def dF(x):
    return -2*x-np.exp(x)

def newton(x, f, df):
    return x - f(x)/df(x)

def iterateNewton(x0, f, df, maxIteration=100, eps=1.5e-16, bound=1000, debug=0):
    x = x0
    xList = [x]
    iterations = maxIteration
    for i in list(range(0,maxIteration)):
        if debug: print("Iteration: ", i+1,"\nx= ", x)
        x = newton(x,f,df)
        xList.append(x)
        if (np.abs(f(x)) < eps or np.abs(x-xList[-2]) < eps):
            iterations = i+1
            break
        if (abs(x-x0) >= bound):
            print("Out of bound! Terminating Newton-Raphson method.")
            iterations = i+1
            break
    return xList, iterations

def Broyden(x0,x1,F,dF,N=100,eps=1.2e-16):
    xList = [x0,x1]
    xn_minus_1 = x0
    xn = x1
    for _ in range(1,N):
        diff_quot =  (xn - xn_minus_1) / (F(xn) - F(xn_minus_1))
        if np.abs(diff_quot) < 1e-16: 
            diff_quot = 1/dF(xn)
            if np.abs(diff_quot) < 1e-16: break
        xn_plus_1 = xn - F(xn) * diff_quot
        xList.append(xn_plus_1)
        if np.abs(F(xn_plus_1)) < eps:
            break
        if np.abs(xn_plus_1 - xn) < eps:
            break
        if np.abs(xn_plus_1 - xn) > 10:
            print('out of bound error')
            break
        xn_minus_1 = xn
        xn = xn_plus_1
    return xList

if __name__ == "__main__":
    x0 = 2.5
    x1 = x0 - F(x0)/dF(x0)

    x_plot = np.arange(-2,2,0.1)
    plt.plot(x_plot,F(x_plot))
    #plt.show()

    # NEWTON
    eps = 1e-16
    error = [F(x0)]
    iterations = 100
    x = x0
    x_newton, iter_newton = iterateNewton(x0,F,dF)
    print('Newton:\nx_root = ', x_newton[-1])
    print(x_newton)
    print(iter_newton)

    x_broyden = Broyden(x0,x1,F,dF)
    print('Broyden')
    print(x_broyden)

    x_star = x_broyden[-1]
    error_broyden = np.abs(np.array(x_broyden[:-1]) - x_star)
    error_newton = np.abs(np.array(x_newton[:-1]) - x_star)

    order_broyden = []
    for xn_plus_1, xn in zip(error_broyden[1:],error_broyden[:-1]):
        if xn_plus_1 > 1e-16 and xn > 1e-16:
            order_broyden.append(np.log(xn_plus_1)/np.log(xn))
    order_newton = []
    for xn_plus_1, xn in zip(error_newton[1:],error_newton[:-1]):
        if xn_plus_1 > 1e-16 and xn > 1e-16:
            order_newton.append(np.log(xn_plus_1)/np.log(xn))
    order_broyden = np.array(order_broyden)
    order_newton = np.array(order_newton)
    print('Order: ',order_broyden)

    plt.clf()
    plt.suptitle('|f(xn)|')
    plt.semilogy(np.arange(0,iter_newton+1),np.abs(F(np.array(x_newton))),"rx-", label='Newton')
    plt.semilogy(np.arange(0,len(x_broyden)),np.abs(F(x_broyden)),"bp-",label='Broyden')
    plt.grid()
    plt.legend()
    scriptpath = os.path.dirname(__file__)
    plt.savefig(os.path.join(scriptpath,"Ex12_1_f(x).png"))

    plt.clf()
    plt.suptitle('Error |x* - xn')
    plt.semilogy(np.arange(0,error_newton.size),error_newton,"rx-", label='Newton')
    plt.semilogy(np.arange(0,error_broyden.size),error_broyden,"bp-",label='Broyden')
    plt.grid()
    plt.legend()
    scriptpath = os.path.dirname(__file__)
    plt.savefig(os.path.join(scriptpath,"Ex12_1_error.png"))  

    plt.clf()
    plt.suptitle('Order of convergence')
    plt.plot(np.arange(2,order_newton.size),order_newton[2:],"rx-", label='Newton')
    plt.plot(np.arange(2,order_broyden.size-1),order_broyden[2:-1],"bp-",label='Broyden')
    plt.grid()
    plt.legend()
    scriptpath = os.path.dirname(__file__)
    plt.savefig(os.path.join(scriptpath,"Ex12_1_order.png")) 

    plt.clf()
    plt.suptitle('Acc as function of evals')
    plt.semilogy(3*np.arange(0,error_newton.size),error_newton,"rx-", label='Newton')
    plt.semilogy(2*np.arange(0,error_broyden.size),error_broyden,"bp-",label='Broyden')
    plt.grid()
    plt.legend()
    scriptpath = os.path.dirname(__file__)
    plt.savefig(os.path.join(scriptpath,"Ex12_1_acc.png"))





