'''
A,B are nxn-matrices
where n = 2^m for simplicites sake

Goal: Compute C = A B
Method: Decompose matrices into 4 blocks of size 2^(m-1)

A =   ( a11 a12)      B = (b11 b12)       C = (c11 c12)  
      ( a21 a22)          (b21 b22)           (c21 c22)
Then the standard multiplication C = A B can be written as follows:
C = (c11 c12) = (a11*b11 + a12*b21      a11*b21 + a12*b22)
    (c21 c22)   (a21*b11 + a22*b21      a21*b21 + a22*b22)
'''

import numpy as np 
import matplotlib.pyplot as plt 
import os 
import time

def matMult(A,B,n):
    C = np.empty_like(A)
    for i in range(n):
        for j in range(n):
            sum = 0.0
            for k in range(n):
                sum = sum + A[i,k]*B[k,j]
            C[i,j] = sum
    return C

def strassen(n,A,B,C,W):
    if n == 1:
        C[0,0] = A*B
        print(C)
        return 0
    n = int(n)
    h = int(n/2)
    print(h,n)
    W[0:h,0:h] = A[0:h,0:h] + A[h:n,h:n]
    print(W[0:h,0:h])
    W[0:h,h:n] = B[0:h,0:h] + B[h:n,h:n]
    print(W[0:h,h:n])
    M1 = np.empty((h,h))
    M1 = strassen(n/2,W[0:h,0:h],W[0:h,h:n],W[h:n,0:h],W[h:n,h:n])
    return W[h:n,0:h]

# def matMultT(A,B,n):
#     C = np.empty_like(A)
#     A = A.transpose()
#     for i in range(n):
#         for j in range(n):
#              sum = 0.0
#             for k in range(n):
#                 sum = sum + A[k,i]*B[k,j]
#             C[i,j] = sum
#     A = A.transpose()
#     return C

if __name__ == "__main__":
    n = 4
    A = np.ones((n,n))
    B = np.ones((n,n)) * 3
    C = np.empty((n,n))
    W = np.empty((n,n))

    start = time.time()
    Cref = np.dot(A,B)
    stopRef = time.time()
    C1 = matMult(A,B,n)
    stop1 = time.time()
    # C2 = matMultT(A,B,n)
    # stop2 = time.time()
    # print("Cref:\n",Cref)
    print("Time REF: ",stopRef - start,"\n","_"*30)
    # print("C1:\n",C1)
    print("matMult close? ",np.allclose(Cref,C1))
    print("Time matMult: ",stop1 - stopRef,"\n","_"*30)


    print("_"*30)
    print("")
    print("Strassen test:")
    print("M1: ",strassen(n,A,B,C,W))
    # print("C2:\n",C2)
    # print("matMultT close? ",np.allclose(Cref,C2))
    # print("Time matMultT: ",stop2 - stop1,"\n","_"*30)