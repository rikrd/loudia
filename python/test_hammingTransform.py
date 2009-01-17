#!/usr/bin/env python

import scipy
import pylab
import ricaudio

a = ricaudio.hammingTransform(50, 1, 1024, 1024, 4)

pylab.plot(a[0,:])

pylab.show()
