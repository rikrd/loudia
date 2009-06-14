#!/usr/bin/env python

import pylab
import scipy

lowBand = 0
highBand = 23
sampleRate = 44100.0
spectrumSize = 2**14

import loudia
m = loudia.BarkBands(lowBand, highBand, sampleRate, spectrumSize)
starts = m.starts()[:,0]

pylab.figure()
for band in range(m.bands()):
    pylab.hold(True)
    
    weight = m.bandWeights( band ).T
    start = starts[band]
    
    x = scipy.arange(start-1, start + weight.shape[1]+1)
    y = [0] + list(weight[0,:])  + [0]

    pylab.plot(x, y, color = 'black')

    ax = pylab.gca()

    # Show half of the spectrum
    ax.set_xlim([0, spectrumSize / 2])
    ax.set_ylim([0.0, 1.1])
    
    # Set the ticks units to radians per second
    ticks = ax.get_xticks()
    ax.set_xticklabels(['%.2f' % (float(tick) / spectrumSize) for tick in ticks])

    # Set the title and labels
    pylab.title('Magnitude of the Frequency Response of a \n Bark Bands implementation')
    pylab.xlabel('Normalized Frequency')
    pylab.ylabel('|H(w)| (no unit)')
    
pylab.show()
