#!/usr/bin/env python

import scipy
import pylab
import loudia

peakSynth = loudia.PeakSynthesize(1024, 4*1024, loudia.Window.HAMMING)

pylab.ion()
for i in range(10):
    a = peakSynth.process(scipy.array([50 + i/10.0, 150 + i/3.0]), scipy.array([-20, -23]))
    pylab.plot(a[0,:])
    
pylab.ioff()

pylab.show()
