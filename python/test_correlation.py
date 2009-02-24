#!/usr/bin/env python

import ricaudio
import scipy

plot = False

atol = 1e-1
rtol = 1e-1

sizeA = 3
a = scipy.reshape(scipy.array(scipy.arange( sizeA ), dtype = 'f4'), (1, sizeA))

sizeB = 10
b = scipy.reshape(scipy.array(scipy.arange( sizeB ), dtype = 'f4'), (1, sizeB))

d = ricaudio.Correlation( sizeA, sizeB )
r = d.process(a, b)
s = scipy.correlate(a[0,:], b[0,:], 'full')
print scipy.allclose(r[0,:], s, atol = atol, rtol = rtol)

d = ricaudio.Correlation( sizeA, sizeB, sizeA + sizeB, -(sizeA + sizeB), True )
r = d.process(b, a)
s = scipy.correlate(a[0,:], b[0,:], 'full')
print scipy.allclose(r[0,:], s, atol = atol, rtol = rtol)

if plot:
    import pylab
    pylab.figure()
    pylab.plot(r[0,:], label = 'ricaudio')
    pylab.plot(s, label = 'scipy')
    pylab.legend()
    pylab.show()

