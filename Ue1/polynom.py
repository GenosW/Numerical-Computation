import numpy as np
from neville import neville

if __name__ == "__main__":
    # Write and format input data (knots and function values)
    data = np.array([[0,0],[1,2],[4,8]], dtype=int)
    x = np.array(data[:,0])
    deg = x.size -1
    y = data[:,1]

    # Polynomial fit
    p_np = np.polyfit(x,y, deg)
    p_np_nice = np.array(p_np)
    # Beautify fit for output
    i = 0
    for n in p_np_nice:
        if abs(n) < 0.001:
            p_np_nice[i] = 0
        i +=1
    p_np = np.poly1d(p_np)
    p_np_nice = np.poly1d(p_np_nice)

    # Compute with neville scheme for comparison
    q, nev = neville(x,y,point=2, retVal = True)
    # q = neville(x,y,point=2)
    # Output
    print('Interpolation polynom for data points:')
    print(data)
    print('x =',x,'\ny =',y)
    print('deg =', deg)
    print('\n-------Interpolation--------')
    print('polyfit: \n',p_np)
    print('\np(x) = ',p_np_nice)
    print('p(2) =', p_np_nice(2))
    print('\n-------Neville scheme-------')
    print(q)
    try: 
        print('nev(2) =',nev)
    except NameError:
        dummy = 0
    print('-------End of program-------')
    