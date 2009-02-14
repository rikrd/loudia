#!/usr/bin/env python

import ricaudio
import scipy

plot = True

npoints = 100
dur = 1.0
f = 4.0
t = scipy.linspace(0, dur, npoints)
x = scipy.cos(2.0 * f * scipy.pi * t)
x = scipy.array(x, dtype = 'f4')
ratio = 1./3.

d = ricaudio.Resample(npoints, npoints/3, ratio, ricaudio.Resample.SINC_BEST_QUALITY)
b = d.process(x)

if plot:
    import pylab
    pylab.subplot(211)
    pylab.plot(x)
    pylab.subplot(212)
    pylab.plot(b[0,:])
    pylab.show()
    
 
#d = ricaudio.Resample(50, 100, 0.5, ricaudio.Resample.SINC_BEST_QUALITY)
#b = d.process(x)

print b
