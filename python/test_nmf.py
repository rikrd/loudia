#!/usr/bin/env python

import ricaudio
import scipy
import pylab

d = ricaudio.NMF(8, 2)
a = scipy.zeros((14, 8), dtype = 'f4')
a[:4, 2] = 1
a[5:9, 6] = 1
a[10:, 2] = 1

print a

components, gains = d.process(a)

print 'components:', components
print 'gains:', gains

pylab.figure()
pylab.plot(components)

pylab.figure()
pylab.plot(gains)
