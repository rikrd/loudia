#!/usr/bin/env python

"""
Levinson-Durbin recursion to solve

Ax = b

where A is a Toelpitz matrix.
"""

import scipy
import scipy.linalg

r = [1, 0.5, 0.5, 0.25]
A = scipy.linalg.toeplitz( r[:-1] )
b = -scipy.array( r[1:] )

def levinson(A, y, ncoeffs = scipy.inf):
    ncoeffs = min(ncoeffs, A.shape[1])

    bs = scipy.zeros( ncoeffs )
    fs = scipy.zeros( ncoeffs )

    bs[-1] = 1 / A[0,0]
    fs[0] = 1 / A[0,0]

    x = scipy.zeros( ncoeffs )
    x[0] = y[0] / A[0,0]
    error = 0

    for i in range(1, ncoeffs):
        ef = scipy.dot( fs[:i] , A[-1, 0:i  ] )
        eb = scipy.dot( bs[-i:], A[0 , 1:i+1] )
        
        den = 1. - (eb * ef)
        
        temp = ( 1. / den * fs[:i+1] - \
                 ef / den * bs[-i-1:]  )
        
        bs[-i-1:] = ( 1. / den * bs[-i-1:] - \
                      eb / den * fs[:i+1] )

        fs[:i+1] = temp
        
        error = scipy.dot( A[-1, :i], x[:i] )
        
        x[:i+1] += ( y[i] - error ) * bs[-i-1:]

    return (x, den)


def levinson2(A, y, ncoeffs = scipy.inf):
    A = scipy.array(A, dtype = 'f4')
    ncoeffs = min(ncoeffs, A.shape[1])
    
    a = scipy.zeros( ncoeffs )
    b = scipy.zeros( ncoeffs )
    
    a[0] = 1
    b[-1] = 1
    
    eps = A[0,0]
    
    x = scipy.zeros( ncoeffs )
    x[0] = y[0] / eps
    
    for i in range(1, ncoeffs):        
        e = - (1. / eps) * scipy.dot(A[-1, :(i)], a[:(i)])
        v = - (1. / eps) * scipy.dot(A[0, 1:(i+1)], b[-(i):])
        
        temp = a[ :(i+1)] + e * b[-(i+1):]
        b[-(i+1):] += v * a[ :(i+1)]
        a[ :(i+1)] = temp
        
        eps *= (1. - (e * v))
        
        ln = y[i] - scipy.dot(A[-1, :(i)], x[:(i)])
        
        x[:(i+1)] += (ln / eps) * b[-(i+1):]
        
    return (x, scipy.sqrt(eps))


def levinsonSymmetricConstrainA1(A, y, ncoeffs = scipy.inf):
    A = scipy.array(A, dtype = 'f4')
    ncoeffs = min(ncoeffs, A.shape[1])
    
    eps = A[0,0]

    ref = scipy.zeros( ncoeffs - 1 )
    
    x = scipy.zeros( ncoeffs )
    x[0] = 1
    
    for i in range(1, ncoeffs):
        gamma = A[0, i] + scipy.dot(A[0, 1:i], x[1:i])

        ref[i-1] = - gamma / eps

        eps *= (1 - ref[i-1]*scipy.conjugate(ref[i-1]))
        
        temp = scipy.conjugate(scipy.flipud(x[:i+1]))
        
        x[:i+1] += ref[i-1] * temp
        x[i] = ref[i-1]
        
    return (x, scipy.sqrt(eps))

    
x, error = levinson2(A, b)

print 'Final error =', error
print 'Final solution x =', x
print 'Final solution A*x =', scipy.dot(A, x)
print 'Expected A*x =', b
print scipy.allclose(scipy.dot(A, x), b)

x, error = levinson(A, b)

print 'Final error =', error
print 'Final solution x =', x
print 'Final solution A*x =', scipy.dot(A, x)
print 'Expected A*x =', b
print scipy.allclose(scipy.dot(A, x), b)

x, error = levinsonSymmetricConstrainA1(A, b)

print 'Final error =', error
print 'Final solution x =', x
print 'Final solution A*x =', scipy.dot(A, x)
print 'Expected A*x =', b
print scipy.allclose(scipy.dot(A, x), b)
