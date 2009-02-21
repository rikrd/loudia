#!/usr/bin/env python

import ricaudio
import scipy

plot = True

a = scipy.array([[1,2,3,4,5,6,7,8,9]], dtype = 'f4')

d = ricaudio.Autocorrelation(9)
r = d.process(a)
s = scipy.correlate(a[0,:], a[0,:], 'full')
s = s[s.shape[0]/2:]
print scipy.allclose(r[0,:], s)



size = 1024
a = scipy.reshape(scipy.array(scipy.arange( size ), dtype = 'f4'), (1, size))
print a.shape

d = ricaudio.Autocorrelation( 1024 )
r = d.process(a)
s = scipy.correlate(a[0,:], a[0,:], 'full')
s = s[s.shape[0]/2:]


print scipy.allclose(r[0,:], s)

if plot:
    import pylab
    pylab.figure()
    pylab.plot(r[0,:], label = 'ricaudio')
    pylab.plot(s, label = 'scipy')
    pylab.legend()
    pylab.show()

