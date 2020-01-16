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

def find_convergence(array):
    index = 0
    for n, np in zip(array[:-1],array[1:]):
        index += 1
        if abs(n - np) < 1e-15:
            break
    return index

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
    l_max = 10000
    index_RQ = []
    index_RQ_with_shift = []
    n_list = []
    for j,l_max in zip(list(range(2,9)),[100,200,600,2100,7500,30000,60000]):
        n = 2**j
        l_max = (100*int(np.power(n,1.2)))

        print('n= {}'.format(n))
        print('l_max= {}'.format(l_max))
        n_list.append(n)
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

        index_RQ.append(find_convergence(error_RQ))
        index_RQ_with_shift.append(find_convergence(error_RQ_with_shift))
        print(index_RQ[-1])
        print(index_RQ_with_shift[-1])


        if j == 3 or j == 5:
            plt.clf()
            plt.suptitle('error for n= {}'.format(n))
            plt.semilogy(np.arange(1,l_max+1), error_RQ, 'r-', label='RQ')
            plt.semilogy(np.arange(1,np.size(error_RQ_with_shift)+1), error_RQ_with_shift, 'b-', label='RQ with shift')
            plt.grid()
            plt.legend()
            plt.savefig(os.path.join(scriptpath,"Ex12_4c_n={}.png".format(n)))

    plt.clf()
    plt.suptitle('number of QR-steps')
    plt.loglog(n_list, index_RQ, 'rx-', label='RQ')
    plt.loglog(n_list, index_RQ_with_shift, 'bx-', label='RQ with shift')
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(scriptpath,"Ex12_4d.png"))