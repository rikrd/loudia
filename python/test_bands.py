#!/usr/bin/env python

lowFreq = 33.0
highFreq = 22000.0
nBands = 34
samplerate = 44100.0
spectralLength = 2048

import ricaudio
m = ricaudio.MelBands(lowFreq, highFreq, nBands, samplerate, spectralLength)
starts = m.starts()

pylab.figure()
for band in range(m.bands()):
    pylab.hold(True)
    weight = m.weight( band )
    print weight
    start = starts[0]
    
    x = scipy.arange(start, start + weight.shape[1])
    y = weight[0,:]

    pylab.plot(x, y)

pylab.show()
