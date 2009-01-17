#!/usr/bin/env python

import scipy
import pylab
import ricaudio

peakSynth = ricaudio.PeakSynthesize(1024, 1024, ricaudio.Window.HAMMING)

pylab.ion()
for i in range(10):
    a = peakSynth.process([50 + i/10.0], [1])
    pylab.plot(a[0,:])
    
pylab.ioff()

pylab.show()
