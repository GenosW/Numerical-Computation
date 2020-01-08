import numpy as np 
import matplotlib.pyplot as plt 
import os 

def R(x,A):
    Ax = np.dot(A,x)
    return np.dot(x.T,Ax) / np.dot(x.T,x)

def g(x,A):
    Ax = np.dot(A,x)
    x2 = np.dot(x.T,x)
    return 2*(x2* Ax + np.dot(x.T,Ax) * x) / x2**2

def getDirection(x,A):
    return -g(x,A)

def armijo(x,d,A,q,sigma):
    t = 1
    for i in range(100):
        xp1 = x + t*d
        if R(xp1,A) < R(x,A) + t*sigma*np.dot(d,d): break
        t *= q
    return t

def checkMin(xn,xnp1):
    return R(xnp1,A) > R(xn,A)

def descent(x0,A,maxSteps):
    q = 0.3
    sigma = 0.99999999999
    eps = 1e-16
    xList = [x0]
    xn = x0
    lambdas = [R(xn,A)]
    steps = maxSteps
    for n in range(maxSteps):
        dn = getDirection(xn,A)
        #print('dir: ',direction)
        stepLength = armijo(xn,dn,A,q,sigma)
        xnp1 = xn + stepLength * dn
        #print(xnp1)
        # if checkMin(xn,xnp1): 
        #     steps = n
        #     break
        lambdas.append(R(xnp1,A))
        xList.append(xnp1)
        xn = xnp1
        # if lambdas[-1] - lambdas[-2] < eps:
        #     steps = n
        #     break
        if n%10 == 0:
            print(n)
            print('lambda = ',lambdas[-1])
            print('stepLength: ',stepLength)
    return xList,lambdas,steps


# n = int(input("n = "))
# print("n = ",n," selected!")
n = 10
h = 1/n
A = 1/h**2 * (2*np.eye(n,n,dtype=float) - np.eye(n,n,k=1,dtype=float) - np.eye(n,n,k=-1,dtype=float))
x0 = np.ones(n)
#x0 = np.random.rand(n)

print(A)
print(x0)
maxSteps = 100
xList,lambdas,steps = descent(x0,A,maxSteps)
print(xList[-1])
EWs, EVs = np.linalg.eigh(A)

print('Numpy: ',EWs[0])
print('Custom: ',lambdas[-1])

# plt.semilogy(range(steps+1),lambdas)
# plt.show()