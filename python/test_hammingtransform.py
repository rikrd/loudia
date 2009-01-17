#!/usr/bin/env python

import scipy
import pylab
import ricaudio

pylab.ion()
for i in range(10):
    a = ricaudio.hammingTransform(50 + i/10.0, 1, 1024, 1024)
    pylab.plot(a[0,:])
    
pylab.ioff()

pylab.show()
