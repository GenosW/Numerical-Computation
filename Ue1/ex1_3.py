import math as m
import matplotlib.pyplot as plt
import numpy as np

from neville import neville

def printLine(length=40, endl='\n'):
    """Prints a line into the console
    """
    line = '-'*length
    # line *= length
    if endl=='\n':
        print(line)
    else:
        print(line, end=endl)

def printStart():
    """Prints the start message
    """
    printLine(length=17, endl='')
    print('START', end='')
    printLine(length=18)

def func(h):
    """ The function we are suppossed to compute\n
    f(h) = (exp(h) - 1)/h
    """
    return (m.exp(h) - 1)/h

def h_i(i):
    """The function to generate the knots\n
    h_i(i) = 2^(-i)
    """
    return 2**(-i)

if __name__ == "__main__":
    # Use neville scheme to approximate f(0) and compute the error
    printStart()
    # Create data points
    # Knots
    I = 10
    iters = list(range(I+1)) # I+1 --> so last element of list is I
    h = list(map(h_i, iters))
    print('Knots:\n',h)
    # Function values
    f = list(map(func,h))
    print('Function values:\n', f)

    # Let's approximate f(0) with the neville scheme
    N, approx = neville(h,f,point=0, retVal=True)
    printLine()
    print('f(0) ≈', approx)
    print('But it should be:\nf(0) = 1')
    printLine()

    # Let's compute the errors
    print('Lets take a look at the errors in a log-log plot:\n*Note that, if error = 1 then it is not part of the scheme anymore')

    fig, ax = plt.subplots()
    b = list(map(lambda x: I+1-x, [0,1,2,3]))
    ax.loglog(h[0:b[0]], abs(N[0:b[0],0]-1), 'x',
            h[0:b[1]], abs(N[0:b[1],1]-1), 's',
            h[0:b[2]], abs(N[0:b[2],2]-1), '^', 
            h[0:b[3]], abs(N[0:b[3],3]-1), 'p')

    ax.set(xlabel='Knots h_i', ylabel='error', title='Error of first 4 neville scheme columns')
    # ax.set(ylim=(10**-6, 1))
    ax.grid()
    # Save fig into file (BEWARE: absolute path!)
    fig.savefig("/Users/peterholzner/Code/NumComp/Ue1/error_nev.png")
    plt.show()
