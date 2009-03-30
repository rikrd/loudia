#!/usr/bin/env python

import loudia
import scipy

plot = False

atol = 1e-1
rtol = 1e-1

size = 9
a = scipy.arange( size )

d = loudia.Autocorrelation( size, size, 0, False )
r = d.process(a)
s = scipy.correlate(a, a, 'full')
s = s[s.shape[0]/2:]
print scipy.allclose(r[0,:], s, atol = atol, rtol = rtol)

d = loudia.Autocorrelation( size, size, 0, True )
r = d.process(a)
s = scipy.correlate(a, a, 'full')
s = s[s.shape[0]/2:]
print scipy.allclose(r[0,:], s, atol = atol, rtol = rtol)

if plot:
    import pylab
    pylab.figure()
    pylab.plot(r[0,:], label = 'loudia')
    pylab.plot(s, label = 'scipy')
    pylab.legend()
    pylab.show()

