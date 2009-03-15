#!/usr/bin/env python

from common import *

plot = True

windowSize = 4
windowHop = 4

d = ricaudio.INMF(8, 2, 10, 0.5, 15, 1e-17)

a = scipy.zeros((14, 8))
a[:4, 2] = 1
a[5:9, 6] = 1
a[10:, 2] = 1

components = []
gains = []
for frame in framer_array(a, windowSize, windowHop):
    c, g = d.process(frame)
    
    gains.append(g)
    components.append(c)

gains = overlap_add(gains, windowSize, windowHop)
components = components[-1]

if plot:
    pylab.figure()
    pylab.imshow(scipy.flipud(a.T))
    pylab.title("Input")
    
    pylab.figure()
    pylab.plot(components.T)
    pylab.title("Components")
    
    pylab.figure()
    pylab.plot(gains)
    pylab.title("Gains")
    
    pylab.show()
