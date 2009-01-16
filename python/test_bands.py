#!/usr/bin/env python

import pylab
import scipy

lowFreq = 133.0
highFreq = 20050.0
nBands = 34
samplerate = 44100.0
spectralLength = 2048

import ricaudio
m = ricaudio.MelBands(lowFreq, highFreq, nBands, samplerate, spectralLength)
starts = m.starts()[:,0]

pylab.figure()
for band in range(m.bands()):
    pylab.hold(True)
    weight = m.bandWeights( band ).T
    print weight.shape
    start = starts[band]
    
    x = scipy.arange(start, start + weight.shape[1])
    y = weight[0,:]

    pylab.plot(x, y)

pylab.show()