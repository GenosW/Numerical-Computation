import math
import numpy
import matplotlib.pyplot as plt
#from IPython import get_ipython

def trapez(interval,N,f):
    
    sum = 0
    h = ( interval[1] - interval[0] ) / N 
    
    for i in range(N):
        sum += ( f(interval[0]+i*h) + f(interval[0]+(i+1)*h) ) * h / 2
        
    return sum


def f1(x):
    
    return x**2


def f2(x):
    
    return abs(x)


def f3(x):
    
    if x < 1/3:
        return 0.5*math.exp(x)
    else:
        return math.exp(x)
    
    
## Code
       
numberPoints = 20
interval = [-1,1]
i = numpy.arange(1,numberPoints+1,dtype=float)
h = 2**-i

error_1 = numpy.zeros(numberPoints)

for i in range(numberPoints):
    N = int((interval[1]-interval[0])/h[i])
    error_1[i] = trapez(interval,N,f1) - 2/3
   
plt.loglog(h,error_1,'bo',label='Error x^2')
plt.xlabel('Intervallbreite')
plt.legend()



error_2 = numpy.zeros(numberPoints)

for i in range(numberPoints):
    N = int((interval[1]-interval[0])/h[i])
    error_2[i] = abs(trapez(interval,N,f2) - 1)
   
plt.plot(h,error_2,'bo',label='Error abs(x)')
plt.xlabel('Intervallbreite')
plt.legend()
plt.show()


error_3 = numpy.zeros(numberPoints)
realValue = 0.5*math.exp(1/3)-0.5*math.exp(-1)+math.exp(1)-math.exp(1/3)

for i in range(numberPoints):
    N = int((interval[1]-interval[0])/h[i])
    error_3[i] =  abs(trapez(interval,N,f3) - realValue)
   
plt.loglog(h,error_3,'ro',label='Error function 3')
plt.xlabel('Intervallbreite')
plt.legend()
plt.show()

