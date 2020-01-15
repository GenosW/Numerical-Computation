import os
import numpy as np 
import matplotlib.pyplot as plt 

def myQR(A,X0,lmax):
    Q, R = np.linalg.qr(X0)
    eigenvalues = [np.diagonal(R)]
    for _ in range(lmax):
        X = np.dot(A,Q)
        Q,R = np.linalg.qr(X)
        eigenvalues.append(np.diagonal(R))
    return X, np.array(eigenvalues)

def myRQ(A,lmax):
    eigenvalues = []
    for _ in range(lmax):
        Q,R = np.linalg.qr(A)
        A = np.dot(R,Q)
        eigenvalues.append(np.diagonal(A))
    return A, np.array(eigenvalues)

def myRQ_with_Shift(A,lmax):
    eps = 1e-14
    eigenvalues = []
    iterations = lmax
    for _ in range(1,lmax+1):
        if np.abs(A[-2,-1]) <= eps * (np.abs(A[-2,-1]) + np.abs(A[-1,-1])):
            shift = 0
        else:
            lambdas = np.real(np.linalg.eigvals(A[-3:-1,-3:-1]))
            shift = lambdas[0]
            if np.abs(lambdas[1] - A[-1,-1]) < np.abs(lambdas[0] - A[-1,-1]):
                shift = lambdas[1]
        I_shift = shift*np.eye(A.shape[0],A.shape[1])
        Q,R = np.linalg.qr(A - I_shift)
        A = np.dot(R,Q) + I_shift
        eigenvalues.append(np.diagonal(A))
        # if np.abs(A[-2,-1]) <= eps * (np.abs(A[-2,-1]) + np.abs(A[-1,-1])):
        #     iterations = i
        #     break
    return A, np.array(eigenvalues), iterations

if __name__ == "__main__":

    scriptpath = os.path.dirname(__file__)
    j = 2
    n = 2**j
    An = 2*np.eye(n,n,k=0) - np.eye(n,n,k=-1) - np.eye(n,n,k=1)
    print(An)
    true_eigenvalues = np.sort(np.linalg.eigvalsh(An))

    print('_'*40)
    print('Basic QR - Algorithm: ')
    X0 = np.eye(n,n,k=0)
    l_max = 10
    A_l_max,eigenvalues_QR = myQR(An,X0,l_max)
    eigenvalues_QR = np.sort(eigenvalues_QR)
    print('A_lmax: ',A_l_max)
    print('Eigenvalues: ',eigenvalues_QR)
    print('True eigenvalues: ', true_eigenvalues)

    print('_'*40)
    print('RQ - Algorithm')
    A_l_max,eigenvalues_RQ = myRQ(An,l_max)
    eigenvalues_RQ = np.sort(eigenvalues_RQ)
    print('A_lmax: ',A_l_max)
    print('Eigenvalues: ',eigenvalues_RQ)
    print('True eigenvalues: ', true_eigenvalues)

    print('_'*40)
    print('RQ - Algorithm with shift')
    A_l_max,eigenvalues_RQ_with_shift, iterations_RQ_with_shift = myRQ_with_Shift(An, l_max)
    eigenvalues_RQ_with_shift = np.sort(eigenvalues_RQ_with_shift)
    print('A_lmax: ',A_l_max)
    print('Eigenvalues: ',eigenvalues_RQ_with_shift)
    print('True eigenvalues: ', true_eigenvalues)
    
    print('_'*50)
    print('PLOTTING')
    l_max = 1000
    for j in [3,5]:
        n = 2**j
        print('n= {}'.format(n))
        An = 2*np.eye(n,n,k=0) - np.eye(n,n,k=-1) - np.eye(n,n,k=1)
        
        true_eigenvalues = np.sort(np.linalg.eigvalsh(An))

        # A_l_max,eigenvalues_QR = myQR(An,X0,l_max)
        # eigenvalues_QR = np.sort(eigenvalues_QR)
        # error_QR = np.max(np.abs(eigenvalues_QR - true_eigenvalues))

        A_l_max,eigenvalues_RQ = myRQ(An,l_max)
        eigenvalues_RQ = np.sort(eigenvalues_RQ)
        error_RQ = np.max(np.abs(eigenvalues_RQ - true_eigenvalues),axis=1)
        #print(error_RQ)

        A_l_max,eigenvalues_RQ_with_shift, iterations_RQ_with_shift = myRQ_with_Shift(An, l_max)
        eigenvalues_RQ_with_shift = np.sort(eigenvalues_RQ_with_shift)
        error_RQ_with_shift = np.max(np.abs(eigenvalues_RQ_with_shift - true_eigenvalues),axis=1)

        plt.clf()
        plt.suptitle('error for n= {}'.format(n))
        plt.semilogy(np.arange(1,l_max+1), error_RQ, 'r-', label='RQ')
        plt.semilogy(np.arange(1,np.size(error_RQ_with_shift)+1), error_RQ_with_shift, 'b-', label='RQ with shift')
        plt.grid()
        plt.legend()
        plt.savefig(os.path.join(scriptpath,"Ex12_4c_n={}.png".format(n)))