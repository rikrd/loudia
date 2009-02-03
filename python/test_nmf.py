#!/usr/bin/env python

import ricaudio
import scipy
import pylab

d = ricaudio.NMF(8, 2)
a = scipy.zeros((10, 8), dtype = 'f4')
components, gains = d.process(a)

print 'components:', components
print 'gains:', gains

pylab.figure()
pylab.plot(components)

pylab.figure()
pylab.plot(gains)
