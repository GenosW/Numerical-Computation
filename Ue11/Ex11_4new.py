import numpy as np 
import matplotlib.pyplot as plt 
import os

def R(x,A):
    val = np.dot(x.T,np.dot(A,x)) / (np.dot(x.T,x))
    return val

def g(x,A):
    x2 = np.dot(x.T,x)
    Ax = np.dot(A,x)
    val =  Ax* x2 + np.dot(np.dot(x.T,Ax),x)
    return 2*val/(x2**2)

def g2(x,A):
    x2 = np.dot(x.T,x)
    Ax = np.dot(A,x)
    val =  x2 * Ax + np.dot(x.T,Ax) * x
    return 2*val/(x2**2)

def armijo(x,d,A,q,sigma):
    t = 1
    for i in range(100):
        t = q**i
        xp1 = x + t*d
        if R(xp1,A) < R(x,A) - t*sigma*np.dot(d,d): break
    return t


fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
axes = np.reshape(axes,4)

for q,ax in zip([0.8,0.5,0.3,0.1],axes):
    for n in [10,20,40,80]:
        h = 1/n
        A = 1/h**2 * (2*np.eye(n,n,dtype=float) - np.eye(n,n,k=1,dtype=float) - np.eye(n,n,k=-1,dtype=float))
        EWs, EVs = np.linalg.eig(A)
        x = np.ones(n)
        maxSteps = 100
        #q = 0.8
        sigma = 1e-13
        lambdas = [R(x,A)]
        for i in range(maxSteps):
            d = -g(x,A)
            tau = armijo(x,d,A,q,sigma)
            x = x + tau * d
            lambdas.append(R(x,A))
        errors = np.array(lambdas) - min(EWs)
        ax.semilogy(range(maxSteps+1),errors,label=str(n))
        print('numpys: ',min(EWs))
        print('custom: ',lambdas[-1])
    ax.grid()
    #ax.legend()
    ax.set_title('q = '+str(q),fontsize=10)
axes[0].set_ylabel("error")
axes[2].set_ylabel("error")
axes[2].set_xlabel("iterations")
axes[3].set_xlabel("iterations")
fig.tight_layout(pad=2, w_pad=0.5, h_pad=0.5)
fig.text(0.04, 0.04, 'sigma =' +str(sigma), fontsize=12)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles,labels,loc='upper left', borderaxespad=0.)
fig.suptitle('Steepest descent',fontsize=14, fontweight='bold')
scriptpath = os.path.dirname(__file__)
plt.savefig(os.path.join(scriptpath,"Ex11_4.png"),dpi=200)

#hab meinen steepest descent doch noch zum laufen gebracht und ein paar verschiedene werte für q durchprobiert