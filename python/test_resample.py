#!/usr/bin/env python

import loudia
import scipy

plot = True

npoints = 100
dur = 1.0
f = 4.0
t = scipy.linspace(0, dur, npoints)
x = scipy.cos(2.0 * f * scipy.pi * t)

ratio = 1./6.
resampleType = loudia.Resample.SINC_BEST_QUALITY

d = loudia.Resample(npoints, int(npoints*ratio), ratio, resampleType)
b = d.process(x)

if plot:
    import pylab
    pylab.figure()
    pylab.subplot(211)
    pylab.plot(x)
    pylab.subplot(212)
    pylab.plot(b[0,:])
    



npoints = 50
dur = 1.0
f = 4.0
t = scipy.linspace(0, dur, npoints)
x = scipy.cos(2.0 * f * scipy.pi * t)
ratio = 3.
resampleType = loudia.Resample.SINC_BEST_QUALITY
 
d = loudia.Resample(npoints, int(npoints*ratio), ratio, resampleType)
b = d.process(x)

if plot:
    import pylab
    pylab.figure()
    pylab.subplot(211)
    pylab.plot(x)
    pylab.subplot(212)
    pylab.plot(b[0,:])
    pylab.show()
