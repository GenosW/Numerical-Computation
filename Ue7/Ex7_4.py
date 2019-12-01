import numpy as np 
import os
import scipy.linalg as la
import matplotlib.pyplot as plt

if __name__ == "__main__":
    scriptpath = os.path.dirname(__file__)
    picpath = os.path.join(scriptpath, "Ex7_4_patterns.pdf")
    n = 10
    # mattype = float
    mattype = int
    e = np.ones(n,dtype=mattype)
    b = np.zeros_like(e)
    b[0] = 1
    I = np.identity(n, dtype=mattype)
    A = 10*I + np.outer(b,e) + np.outer(e,b)
    A2 = A[n::-1,n::-1] 

    # print(np.outer(b,e))
    # print(np.outer(e,b))
    # print(b)
    # # print(I)
    print(A)
    # print(A2)
    P, L, U = la.lu(A)
    P2, L2, U2 = la.lu(A2)
    fig, axs = plt.subplots(2,3)
    ax1 = axs[0,0]
    ax2 = axs[1,0]
    ax3 = axs[0, 1]
    ax4 = axs[1, 1]
    ax5 = axs[0, 2]
    ax6 = axs[1, 2]
    ax1.spy(A)
    ax2.spy(A2)
    ax3.spy(L)
    ax4.spy(L2)
    ax5.spy(U)
    ax6.spy(U2)
    fig.savefig(picpath)
    # plt.show()