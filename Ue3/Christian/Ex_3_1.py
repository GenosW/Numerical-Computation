import numpy as np
import math
import matplotlib.pyplot as plt

def cheby(interval,n):
    
    a = interval[0]
    b = interval[1]
    i = np.arange(0,n+1)
    xicheb = np.cos(math.pi*(2*i+1)/(2*n+2))
    xi = (a+b)/2 + (b-a)/2*xicheb
    
    return xi


    
def neville(x,f,val):

    M = np.zeros([x.size,x.size])
    M[:,0] = f
    
    for m in range(1,x.size):
        for j in range(0,x.size-m):
            M[j,m] = ( (val-x[j]) * M[j+1,m-1] - (val-x[j+m]) * M[j,m-1] ) / (x[j+m] - x[j] )
            

    return M[0,x.size-1]


def Taylor(x,n):
    
    sum = 0
    for i in range(0,n+1):
        sum = sum + np.power(x,(2*i)) / np.power(4,(i+1))
        
    return sum


def f(x):
    return 1/(4-np.power(x,2))

## Code
    
errorpoints = 100
interval = [-1, 1]
degreemax = 20

error = np.zeros(degreemax)
errorTaylor = np.zeros(degreemax)
c = 0

for n in np.arange(1,degreemax+1):
    
    
    errortemp = np.zeros(errorpoints)
    xi = cheby(interval,n)
    fi = f(xi)
    x = np.linspace(-1,1,errorpoints)
    
    for i in range(errorpoints):
        errortemp[i] = neville(xi,fi,x[i])-f(x[i])
        
    error[c] = abs(max(errortemp))
    
    if(n%2 == 0):
        errorTaylor[c] = abs(Taylor(1,n) - f(1))
    else:
        errorTaylor[c] = abs(Taylor(1,n-1) - f(1))
        
    c += 1

plt.semilogy(np.arange(1,degreemax+1),error,'bo', label='Error Cheby') 
plt.semilogy(np.arange(1,degreemax+1),errorTaylor,'ro', label='Error Taylor')
plt.legend() 
