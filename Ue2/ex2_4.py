import math as m
import time

import numpy as np


def S(N, mode = 'small', a = [0,0,0]):
    if N>500 and mode == 'small':
        mode = 'medium'
    if mode == 'small':
        if N > 0:
            return 1/N**2 + S(N-1)
        else:
            return N
    elif mode == 'medium':
        i = 1
        sum = 1
        while(i and i<N):
            i += 1
            sum += 1/i**2
        return sum
    elif mode == 'large':
        val = a[0] + a[1]/N + a[2]/(N**2)
        return val
    else:
        print('Not a mode!')
        return 0

def normalize(*args):
    arglist = []
    for i in args: arglist.extend(i)
    minimum = min(arglist)
    for i in range(len(arglist)): arglist[i] /= minimum
    return arglist, minimum

if __name__ == "__main__":
    S_small = []
    S_large = []
    t_small = []
    t_large = []
    errors_large = []
    errors_small = []
    B = []
    C = []
    # C*a=B
    N = (10**2, 10**3, 10**4, 10**5)
    N_numbers = len(N)
    for z in N:
        start = time.time() 
        S_small.append(S(z))
        t_small.append(time.time()-start)
    # Compute coefficients for large mode
    for z in [10, 100, 1000]: B.append(S(z)*(z**2))
    for j in range(1,4): C.append([100**j, 10**j, 1])
    a = np.linalg.solve(C, B)
    print('System of equations:')
    for j,i,vec in zip(C[:],B,['a0','a1','a2']): print(j,vec,'|',i)
    print('Coefficients: ', a)
    print('Coeff solution correct? :', np.allclose(np.dot(C, a), B))
    # Use approximation
    for z in N: 
        start = time.time() 
        S_large.append(S(z, mode = 'large', a = a))
        t_large.append(time.time()-start)
    # Normalize runtimes to smallest runtime
    arglist, unit = normalize(t_small, t_large)
    # Compute errors
    for i in S_large: errors_large.append(abs(i-(m.pi**2)/6))
    for i in S_small: errors_small.append(abs(i-(m.pi**2)/6))

    # Write results into console
    print('-----------------------------')
    print('For the numbers N=', N)
    print('computing S(N) = ∑(1/n**2) from 1 to N.')
    print('lim S(N for N->inf = (pi^2)/6 =', (m.pi**2)/6)
    print('---------Method1---------')
    print('Straight forward approach')
    print('Results: ', S_small)
    print('Error:   ', errors_large)
    print('Runtimes:', arglist[:N_numbers])
    print('---------Method2---------')
    print('Approximation S(N) ≈ a0 + a1/N + a2/N^2')
    print('Results: ', S_large)
    print('Error:   ', errors_large)
    print('Runtimes:', arglist[N_numbers:])
    print('-------------------------')
    print('runtimes in multiples of',unit)

