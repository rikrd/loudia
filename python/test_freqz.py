#!/usr/bin/env python

import loudia
import scipy

# Create the omega array
npoints = 1000
w = scipy.arange(-scipy.pi, scipy.pi, 2.0*scipy.pi/npoints)

# Create the coefficients of the filter
nbcoeffs = 10
b = scipy.ones((nbcoeffs, 1)) / nbcoeffs
a = scipy.ones((1, 1))

print b
print a
print w

# Calculate the frequency response
d = loudia.freqz(b, a, w)

# Plot the frequency response
import pylab

pylab.subplot(2,1,1)
pylab.plot(w, abs(d[:,0]))
pylab.title('Magnitude of the Frequency Response')

pylab.subplot(2,1,2)
pylab.plot(w, scipy.angle(d[:,0]))
pylab.title('Angle of the Frequency Response')

pylab.show()
