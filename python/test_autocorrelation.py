#!/usr/bin/env python

import ricaudio
import scipy

plot = False

atol = 1e-1
rtol = 1e-1

size = 9
a = scipy.reshape(scipy.array(scipy.arange( size ), dtype = 'f4'), (1, size))

d = ricaudio.Autocorrelation( size, size, 0, False )
r = d.process(a)
s = scipy.correlate(a[0,:], a[0,:], 'full')
s = s[s.shape[0]/2:]
print scipy.allclose(r[0,:], s, atol = atol, rtol = rtol)



size = 128
a = scipy.reshape(scipy.array(scipy.arange( size ), dtype = 'f4'), (1, size))

d = ricaudio.Autocorrelation( size, size, 0, True )
r = d.process(a)
s = scipy.correlate(a[0,:], a[0,:], 'full')
s = s[s.shape[0]/2:]
print scipy.allclose(r[0,:], s, atol = atol, rtol = rtol)

if plot:
    import pylab
    pylab.figure()
    pylab.plot(r[0,:], label = 'ricaudio')
    pylab.plot(s, label = 'scipy')
    pylab.legend()
    pylab.show()

