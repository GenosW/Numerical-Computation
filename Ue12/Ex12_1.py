import numpy as np
import matplotlib.pyplot as plt 

def F(x):
    return 2-x-np.exp(x) 

def dF(x):
    return -1-np.exp(x)

def netwonXD(xn,f,df):
    if np.linalg.det(df(xn)):
        res = f(xn)
        corr = np.linalg.solve(df(xn),res)
        return xn - corr, res
    else:
        print("df is noninvertible->cant be solved...returning input point x_n!")
        return xn

def Broyden(x0,x1,N=100,eps=1e-16):
    xList = [x0,x1]
    xn_minus_1 = x0
    xn = x1
    for _ in range(1,N):
        denominator = F(xn) - F(xn_minus_1)

        xn_plus_1 = -(xn - xn_minus_1)*F(xn)
        if np.abs(denominator) < 1e-15: 
            xn_plus_1 /= (denominator)
        else:
            xn_plus_1 /= dF(xn)
        xn_plus_1 += xn

        xList.append(xn_plus_1)

        if np.abs(F(xn_plus_1)) <= eps:
            break
    return xList

if __name__ == "__main__":
    x0 = 2.8
    x1 = x0 - F(x0)/dF(x0)

    # NEWTON
    for i in list(range(0,100)):
        xn1, res = netwonXD(x,f,df)
        error.append(np.linalg.norm(res))
        if np.allclose(xn1,x,rtol=eps):
            # print("x:\n",x)
            # print("xn1:\n",xn1)
            print("Break condition: <xn, xnp1 close>")
            iterations = i+1
            np.copyto(x, xn1)
            break
        np.copyto(x, xn1)



