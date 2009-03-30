#!/usr/bin/env python

import loudia
import scipy

plot = False

atol = 1e-1
rtol = 1e-1

sizeA = 3
a = scipy.arange( sizeA )

sizeB = 10
b = scipy.arange( sizeB )

d = loudia.Correlation( sizeA, sizeB )
r = d.process(a, b)
s = scipy.correlate(a, b, 'full')
print scipy.allclose(r[0,:], s, atol = atol, rtol = rtol)

d = loudia.Correlation( sizeA, sizeB, sizeA + sizeB, -(sizeA + sizeB), True )
r = d.process(b, a)
s = scipy.correlate(a, b, 'full')
print scipy.allclose(r[0,:], s, atol = atol, rtol = rtol)

if plot:
    import pylab
    pylab.figure()
    pylab.plot(r[0,:], label = 'loudia')
    pylab.plot(s, label = 'scipy')
    pylab.legend()
    pylab.show()

