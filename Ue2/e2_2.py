import numpy as np
import math as m
import matplotlib.pyplot as plt

from neville import neville

def f(x ):
    return m.tan(x)

def g(x):
    return (max(x,0))**1.5

def D(h, x0 = 0, func = f):
    return (func(x0 + h) - func(x0 - h))/(2*h)

def plotCustom(nev, h, fun, name=''):
    fig, ax = plt.subplots()
    labels = ['col1', 'col2', 'col3']
    if fun == f:
        exact = 1
        labelsh = ['h^2', 'h^4', 'h^8']
        h2 = list(map(lambda x: x**2, h))
        h3 = list(map(lambda x: x**4, h))
        h4 = list(map(lambda x: x**8, h))
    else:
        exact = 0 #0.5**1.5
        labelsh = ['h', 'h^(1/2)', 'h^(1/4)']
        h2 = list(map(lambda x: x, h))
        h3 = list(map(lambda x: m.sqrt(x), h))
        h4 = list(map(lambda x: m.sqrt(m.sqrt(x)), h))
    ax.loglog(h[:], abs(nev[:,0]-exact), 'x', label=labels[0])
    ax.loglog(h[:-1], abs(nev[:-1,1]-exact), 's', label=labels[1])
    ax.loglog(h[:-2], abs(nev[:-2,2]-exact), '^', label=labels[2])
    ax.loglog(h, h2, '-', label=labelsh[0])
    ax.loglog(h, h3, '-', label=labelsh[1])
    ax.loglog(h, h4, '-', label=labelsh[2])
    ax.legend()
    ax.set(xlabel='Knots h_i', ylabel='error', title='Error of first 3 neville scheme columns')
    i = ax.get_ylim()
    if i[0] < 10**(-17):
        ax.set_ylim(bottom=10**(-17))
        print(i)
    # ax.set(ylim=(10**-6, 1))
    ax.grid()
    # Save fig into file (BEWARE: absolute path!)
    fig.savefig("/Users/peterholzner/Code/NumComp/Ue2/"+name)
    # plt.show()

if __name__ == "__main__":
    print("Ex2_2:_______")
    N = 14
    h = list(map(lambda i: 2**(-i), range(N)))
    h2 = list(map(lambda x: x**2, h))
    for fun in [f, g]:
        Dsym = [D(x, func=fun) for x in h]
        Dsym2 = [D(x, func=fun) for x in h2]
        # print(Dsym)
        nev, res = neville(h, Dsym, retVal=True)
        nev2 = neville(h2, Dsym2)
        print('Result =', res)
        print('Matrix:')
        print(nev)
        #Plot Dsym
        if fun == f:
            name1 = 'ex2_2-fh.png'
            name2 = 'ex2_2-fh2.png'
        else:
            name1 = 'ex2_2-gh.png'
            name2 = 'ex2_2-gh2.png'
        plotCustom(nev, h, fun, name=name1)
        plotCustom(nev2, h, fun, name=name2)
        
        # # Plot Dsym2
        # fig, ax = plt.subplots() 
        # ax.loglog(h[:], abs(nev2[:,0]-exact), 'x',
        #         h[:-1], abs(nev2[:-1,1]-exact), 's',
        #         h[:-2], abs(nev2[:-2,2]-exact), '^',
        #         h, h2, '-',
        #         h, h3, '-',
        #         h, h4, '-')
        # ax.set(xlabel='Knots h_i', ylabel='error', title='Error of first 3 neville scheme columns')
        # # ax.set(ylim=(10**-6, 1))
        # ax.grid()
        # # Save fig into file (BEWARE: absolute path!)
        # if fun == f:
        #     fig.savefig("/Users/peterholzner/Code/NumComp/Ue2/ex2_2-fh2.png")
        # else:
        #     fig.savefig("/Users/peterholzner/Code/NumComp/Ue2/ex2_2-gh2.png")
        # # plt.show()
