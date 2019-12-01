import numpy as np
import math


def cholesky(A):
    
    C = np.zeros([len(A),len(A)])
    n = len(A)
    
    for i in range(0,n):
        for k in range(i,n):
            if(k==i):
                sum = 0
                for j in range(0,i):
                    sum += np.power(C[i][j],2)
                C[i][i] = math.sqrt(A[i][i]-sum)
            else:
                sum = 0
                for j in range(0,i):
                    sum += C[i][j]*C[i+1][j]
                C[k][i] = (A[i][k]-sum)/C[i][i]
    
    
    return C


## Code
    
# define an example Matrix A
    
# A = [[1,-1,2],[-1,5,-4],[2,-4,6]]
N = 4
A_base = np.random.random_integers(-10,10,size=(N,N))
A = np.array((A_base + A_base.T)/2)
C = cholesky(A)
C = np.array(C)
A = np.array(A)
print('A is: \n')
print(A)
print('\nC is: \n')
print(C)
print(" ")
print(np.dot(C,C.transpose()))
# print(np.dot(C.transpose(),C))
