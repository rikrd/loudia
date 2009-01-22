#!/usr/bin/env python

import ricaudio
import scipy

npoints = 1000
w = scipy.arange(-scipy.pi, scipy.pi, 2*scipy.pi/(npoints), dtype = 'f4')

nbcoeffs = 10
b = scipy.ones((nbcoeffs, 1), dtype = 'f4') / nbcoeffs

a = scipy.ones((1, 1), dtype = 'f4')

d = ricaudio.freqz(b, a, w)

import pylab
pylab.subplot(2,1,1)
pylab.plot(w, abs(d[:,0]))
pylab.title('Magnitude of the Frequency Response')

pylab.subplot(2,1,2)
pylab.plot(w, scipy.angle(d[:,0]))
pylab.title('Angle of the Frequency Response')

pylab.show()
