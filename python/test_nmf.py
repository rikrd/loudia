#!/usr/bin/env python

import ricaudio
import scipy
import pylab

plot = True

d = ricaudio.NMF(8, 2, 100, 1, 1e-7)
a = scipy.zeros((140, 8), dtype = 'f4')
a[:40, 2] = 1
a[50:90, 6] = 1
a[100:, 2] = 1

print a

components, gains = d.process(a)

print 'components:', components
print 'gains:', gains

if plot:
    pylab.figure()
    pylab.imshow(scipy.flipud(a.T))
    pylab.title("Input")
    
    pylab.figure()
    pylab.plot(components)
    pylab.title("Components")
    
    pylab.figure()
    pylab.plot(gains.T)
    pylab.title("Gains")
    
    pylab.show()
