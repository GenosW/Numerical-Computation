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

M1 = (A11 + A22) (B11 + B22)
M2 = (A21 + A22) B11
M3 = A11 (B12 − B22)
M4 = A22 (B21 − B11)
M5 = (A11 + A12) B22
M6 = (A21 − A11) (B11 + B12)
M7 = (A12 − A22) (B21 + B22)
C11 = M1 + M4 − M5 + M7
C12 = M3 + M5
C21 = M2 + M4
C22 = M1 − M2 + M3 + M6
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
        #print(C)
        return 0
    n = int(n)
    h = int(n/2)
    # print('(n,h) = ({},{})'.format(n,h))
    # Partition... for better readability
    # A
    A11 = A[0:h,0:h]
    A12 = A[0:h,h:n]
    A21 = A[h:n,0:h]
    A22 = A[h:n,h:n]
    # B
    B11 = B[0:h,0:h]
    B12 = B[0:h,h:n]
    B21 = B[h:n,0:h]
    B22 = B[h:n,h:n]
    #C
    C11 = C[0:h,0:h]
    C12 = C[0:h,h:n]
    C21 = C[h:n,0:h]
    C22 = C[h:n,h:n]
    # W
    W11 = W[0:h,0:h]
    W12 = W[0:h,h:n]
    W21 = W[h:n,0:h]
    W22 = W[h:n,h:n]
    
    ## Compute
    # M1 = (A11 + A22) (B11 + B22)
    W11 = A11 + A22
    W12 = B11 + B22
    strassen(h,W11,W12,W22,W21)
    C11 = W22
    C22 = W22

    # M2 = (A21 + A22) B11
    W11 = A21 + A22
    W12 = B11
    strassen(h,W11,W12,W22,W21)
    C21 = W22
    C22 -= W22
    
    # M3 = A11 (B12 − B22)
    W11 = A11
    W12 = (B12 - B22)
    strassen(h,W11,W12,W22,W21)
    C12 = W22
    C22 += W22

    # M4 = A22 (B21 − B11)
    W11 = A22
    W12 = (B21 - B11)
    strassen(h,W11,W12,W22,W21)
    C11 += W22
    C21 += W22

    # M5 = (A11 + A12) B22
    W11 = (A11 + A12)
    W12 = B22
    strassen(h,W11,W12,W22,W21)
    C11 -= W22
    C12 += W22

    # M6 = (A21 − A11) (B11 + B12)
    W11 = (A21 - A11)
    W12 = (B11 + B12)
    strassen(h,W11,W12,W22,W21)
    C22 += W22

    # M7 = (A12 − A22) (B21 + B22)
    W11 = (A12 - A22)
    W12 = (B21 + B22)
    strassen(h,W11,W12,W22,W21)
    C11 += W22

    # print("M7 done -> C:\n",C)
    # if True: return 0
    return 0

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
    n = 2**8
    print('n= ',n)
    A = np.ones((n,n))
    B = np.ones((n,n)) * 3
    C = -10*np.ones((n,n))
    W = np.empty((n,n))

    start = time.time()
    Cref = np.dot(A,B)
    stopRef = time.time()
    C1 = matMult(A,B,n)
    stop1 = time.time()
    startStr = time.time()
    ret = strassen(n,A,B,C,W)
    stopStr = time.time()
    # C2 = matMultT(A,B,n)
    # stop2 = time.time()
    # print("Cref:\n",Cref)
    print("Time numpy: ",stopRef - start)
    #print("Cref:\n",Cref)
    print("_"*30)
    print("Time matMult: ",stop1 - stopRef)
    #print("C1:\n",C1)
    print("matMult close? ",np.allclose(Cref,C1))
    print("_"*30)
    print("Strassen test:")
    print("Strassen: ",ret)
    print("Time Strassen: ",stopStr - startStr)
    #print("C:\n",C)
    print("Strassen close? ",np.allclose(Cref,C))
    # print("C2:\n",C2)
    # print("matMultT close? ",np.allclose(Cref,C2))
    # print("Time matMultT: ",stop2 - stop1,"\n","_"*30)