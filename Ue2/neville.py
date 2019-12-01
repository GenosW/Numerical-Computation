import numpy as np

def neville(x, f, point=0, n=-1, retVal = False):
    """ Neville scheme:\n
    x...        vector of knots of length n+1,\n
    f...        vector of data values of length n+1,\n
    point...    the point(knot) to evaluate the neville scheme at,\n
    n...        size of the neville scheme, calculated automatically if left at default,\n
    retVal...   determines what is given as return value (scheme or scheme + value(point)).
    """
    
    if n<0:
        x = np.array(x)
        n = x.size
    q = np.zeros((n, n))
    q[:,0] = np.copy(f)
    # q[:,0] = f

    for m in range(1, n):
        for i in range(0, n-m):
            q[i,m] = (point - x[i])*q[i+1,m-1] - (point - x[i+m])*q[i,m-1]
            q[i,m] /= (x[i+m] - x[i])
    result = q
    if retVal:
        result = [q,q[i,m]]
    return result

def neville_vec(x, f, point=0, n=-1, retVal = False):
    """ Neville scheme:\n
    x...        vector of knots of length n+1,\n
    f...        vector of data values of length n+1,\n
    point...    the point(knot) to evaluate the neville scheme at,\n
    n...        size of the neville scheme, calculated automatically if left at default,\n
    retVal...   determines what is given as return value (scheme or scheme + value(point)).
    """
    
    if n<0:
        x = np.array(x)
        n = x.size
    q = np.zeros((n, n))
    q[:,0] = np.copy(f)
    # q[:,0] = f

    for m in range(1, n):
        for i in range(0, n-m):
            q[i,m] = (point - x[i])*q[i+1,m-1] - (point - x[i+m])*q[i,m-1]
            q[i,m] /= (x[i+m] - x[i])
    result = q
    if retVal:
        result = [q,q[i,m]]
    return result

if __name__ == "__main__":
    # TEST
    x = np.array([0,1,3])
    f = np.array([1,3,2])
    q, res = neville(x, f, point=2, retVal=True)
    print('f(2) =', res)
    print('Scheme:\n', q)
