<<<<<<< HEAD
import numpy as np 
import matplotlib.pyplot as plt 
import os 

def newtonContinuation(xn,s,H,dH):
    if dH(xn,s):
        res = H(xn,s)
        corr = res/dH(xn,s)
        xNew = xn - corr
        #print('xNew= ',xNew)
        return xNew
    else:
        print('Derivative is zero, aborting newton method and returning last x!')
        return xn

def f(x):
    return np.arctan(x)

def df(x):
    return 1/(1+x*x)

def H(x,s):
    return f(x) - (1-s)*np.arctan(4)

def dH(x,s):
    return df(x)

if __name__ == "__main__":
    print('NumComp Ex 11, Task 1:')
    Si = np.linspace(0,1,num=11)
    print("s_i: ",Si)
    x = 4
    print('Starting guess x_0: ',x)
    print('H(x_0,0)= ',H(x,0))
    print('Getting to work, boss!')
    print('-'*30)

    xi = [x]
    fi = [np.abs(f(x))]
    eps = 1e-16
    iterations = 0
    for s in Si[1:]:
        its = 0
        for i in range(10):
            x = newtonContinuation(x,s,H,dH)
            its += 1 
            if np.allclose(x,xi[-1]) or np.abs(H(x,s))<eps: break
            if s == Si[-1]:
                xi.append(x)
                fi.append(np.abs(f(x)))
        iterations += its
        xi.append(x)
        fi.append(np.abs(f(x))) 
        print('Intermediate result after {} iterations (for this s)'.format(its))
        print('s= {}, x= {}, f(x)= {}'.format(s,x,f(x)))
    
    print("Arrived at solution after ",iterations," iterations.")
    print("final:\n",x)
    #print("final error:\n",np.linalg.norm(np.dot(A,x) - g(x)))
    #x2 = np.power(1e-2,np.arange(0,iterations+1))
    plt.semilogy(xi,fi,"bx-",label="Newton")
    # plt.semilogy(np.arange(0,iterations+1),x2,"b-", label="order x^2 ???")
    plt.grid()
    plt.legend()
    # plt.show()
    scriptpath = os.path.dirname(__file__)
    plt.savefig(os.path.join(scriptpath,"Ex10_3.png"))

=======
import numpy as np 
import matplotlib.pyplot as plt 
import os 

def newtonContinuation(xn,s,H,dH):
    if dH(xn,s):
        res = H(xn,s)
        corr = res/dH(xn,s)
        xNew = xn - corr
        #print('xNew= ',xNew)
        return xNew
    else:
        print('Derivative is zero, aborting newton method and returning last x!')
        return xn

def f(x):
    return np.arctan(x)

def df(x):
    return 1/(1+x*x)

def H(x,s):
    return f(x) - (1-s)*np.arctan(4)

def dH(x,s):
    return df(x)

if __name__ == "__main__":
    print('NumComp Ex 11, Task 1:')
    Si = np.linspace(0,1,num=11)
    print("s_i: ",Si)
    x = 4
    print('Starting guess x_0: ',x)
    print('H(x_0,0)= ',H(x,0))
    print('Getting to work, boss!')
    print('-'*30)

    xi = [x]
    fi = [np.abs(f(x))]
    eps = 1e-16
    iterations = 0
    for s in Si[1:]:
        its = 0
        for i in range(10):
            x = newtonContinuation(x,s,H,dH)
            its += 1 
            if np.allclose(x,xi[-1]) or np.abs(H(x,s))<eps: break
            if s == Si[-1]:
                xi.append(x)
                fi.append(np.abs(f(x)))
        iterations += its
        xi.append(x)
        fi.append(np.abs(f(x))) 
        print('Intermediate result after {} iterations (for this s)'.format(its))
        print('s= {}, x= {}, f(x)= {}'.format(s,x,f(x)))
    
    print("Arrived at solution after ",iterations," iterations.")
    print("final:\n",x)
    #print("final error:\n",np.linalg.norm(np.dot(A,x) - g(x)))
    #x2 = np.power(1e-2,np.arange(0,iterations+1))
    plt.semilogy(xi,fi,"bx-",label="Newton")
    # plt.semilogy(np.arange(0,iterations+1),x2,"b-", label="order x^2 ???")
    plt.grid()
    plt.legend()
    # plt.show()
    
    scriptpath = os.path.dirname(__file__)
    plt.savefig(os.path.join(scriptpath,"Ex11_1.png"))

>>>>>>> 8e80f12d195b7a7b603ea827258841e148ebe0e7
    