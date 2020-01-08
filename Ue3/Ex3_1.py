import numpy as np
import math as m
import matplotlib.pyplot as plt

def f(x):
    return 1/(4-x**2)

def chebyshev(a, b, n):
    i = np.arange(0,n+1)
    xCheb = np.cos(m.pi *(i+0.5)/(n+1))
    return 0.5*(a+b+(b-a)*xCheb)

def neville(x, f, point=0, n=-1, retVal = "point"):
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
    else:
        result = q
    return result

def taylor_of_f(h, n, a=0):
    sum = 0
    for i in np.arange(0,n+1,step=2):
        sum += (h-a)**i / 2**(i+2)
    return sum

if __name__ == "__main__":
    a = -1
    b = 1
    errorPoints = 100
    degreeMax = 10
    x = np.linspace(a,b,errorPoints)
    errorCheby = np.zeros(degreeMax)
    errorTaylor = np.zeros(degreeMax)
    degrees = np.arange(1,degreeMax+1)

    for n in degrees:
        xi = chebyshev(a,b,n)
        fi = f(xi)
        errorTemp = np.zeros(errorPoints)

        for i in range(errorPoints):
            errorTemp[i] = abs(neville(xi,fi,x[i])-f(x[i]))

        errorCheby[n-1] = max(errorTemp)

        for i in range(errorPoints):
            if(n%2 == 0):
                errorTemp[i] = abs(taylor_of_f(x[i],n) - f(x[i]))
            else:
                errorTemp[i] = abs(taylor_of_f(x[i],n-1) - f(x[i]))
        errorTaylor[n-1] = max(errorTemp)

plt.semilogy(np.arange(1,degreeMax+1),errorCheby,'bo', label='Error Cheby') 
plt.semilogy(np.arange(1,degreeMax+1),errorTaylor,'ro', label='Error Taylor')
plt.legend()
plt.grid()
plt.savefig("/Users/peterholzner/Code/NumComp/Ue3/Ex3_1") 
plt.show()
