import numpy as np
import matplotlib.pyplot as plt

def steepestDescentStep(A,x,b):
    r = b - np.dot(A,x)
    alpha = np.dot(r,r)/np.dot(r,np.dot(A,r))
    x = x+alpha*r
    return x


n = 10
h = 1/n
n = n*n
l = 200
A = 2*np.eye(n,n,k=0) - np.eye(n,n,k=-1) - np.eye(n,n,k=1)
A = A / (h*h)
b = np.ones(n)
x_star = np.linalg.solve(A,b)

## Steepest descent
x = np.zeros(n)
errorSD = []
residualSD = []

for i in range(l):
    x = steepestDescentStep(A,x,b)
    residualSD.append(np.linalg.norm(b-np.dot(A,x)))
    errorSD.append(np.linalg.norm(x_star-x))

## Conjugate Gradient (Hestenes & Stiefel)

x = np.zeros(n)
r = b
p = r
errorCG = []
residualCG = []
r_before = np.ones(n)

for i in range(l):
    alpha = np.dot(r,r)/np.dot(p,np.dot(A,p))
    x = x + alpha*p
    np.copyto(r_before,r)
    r = r - alpha*np.dot(A,p)
    beta = np.dot(r,r)/np.dot(r_before,r_before)
    p = r + beta*p
    residualCG.append(np.linalg.norm(b-np.dot(A,x)))
    errorCG.append(np.linalg.norm(x_star-x))


## Plotting the results

plt.semilogy(range(l),errorSD,label='errorSD')

plt.semilogy(range(l),residualSD,label='residualSD')
plt.semilogy(range(l),errorCG,label='errorCG')
plt.semilogy(range(l),residualCG,label='residualCG')
plt.grid()
plt.legend()
plt.show()
