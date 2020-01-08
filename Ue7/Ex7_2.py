import numpy as np

def cholesky(A):
    n, m = A.shape
    C = np.zeros_like(A, dtype=float)
    if (m!=n):
        return C
    for i in np.arange(n-1, -1, -1):
        if i==n-1:
            C[i,i] = np.sign(A[i,i])*np.sqrt(abs(A[i,i]))
        # elif i==n-2:
        #      #- sum([np.power(C[k,i],2) for k in np.arange(i+1,n,1)]))
        #     C[i,i] =  np.sqrt(A[i,i] - C[i+1,i]**2)
        else:
            summe = 0.0
            for k in np.arange(i+1,n):
                summe += C[k,i]**2
            C[i,i] =  np.sign(A[i,i])*np.sqrt(abs(A[i,i] - summe))
        for j in np.arange(i-1, -1, -1):
            C[i,j] = (A[i,j] - sum([C[k,i]*C[k,j] for k in np.arange(max((i,j)),n,1)]))/C[i,i]
    return C

if __name__ == "__main__":
    A = np.array([[1,-1,2],[-1,5,-4],[2,-4,6]])
    # N = 4
    # A_base = np.random.random_integers(-10,10,size=(N,N))
    # A = np.array((A_base + A_base.T)/2)
    C = cholesky(A)

    print('A is: \n')
    print(A)
    print('\nC is: \n')
    print(C)
    print("Ct*C is:")
    print(np.dot(C.transpose(), C))