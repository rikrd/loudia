#!/usr/bin/env python

import ricaudio
from sepel.inputs import pyricaudio
import pylab
import sys
import scipy

filename = sys.argv[1]

frameSize = 1024 / 44100.0
frameStep = 512 / 44100.0

fftSize = 1024

analysisLimit = 10.0

# Creation of the pipeline        
stream = pyricaudio.sndfilereader({'filename': filename,
                                   'windowSizeInTime': frameSize,
                                   'windowStepInTime': frameStep,
                                   'encodingKey': 'encoding',
                                   'channelCountKey': 'channelCount',
                                   'samplesOriginalKey': 'samples',
                                   'samplesKey': 'samplesMono',
                                   'samplerateKey': 'samplerate',
                                   'timestampBeginKey': 'timestampBegin',
                                   'timestampEndKey': 'timestampEnd',
                                   'limit':analysisLimit})


stream = pyricaudio.window_ricaudio(stream, {'inputKey': 'samplesMono',
                                             'outputKey': 'windowed',
                                             'windowType': 'blackmanharris'})

stream = pyricaudio.fft_ricaudio(stream, {'inputKey': 'windowed',
                                          'outputKey': 'fft',
                                          'zeroPhase': True,
                                          'fftLength': fftSize})


#plots = ['mag', 'phase', 'peaks']

plots = ['mag', 'phase', 'peaks']

plotSize = fftSize / 4

pylab.ion()

if 'peaks' in plots:
    peaker = ricaudio.PeakPick( 10 )

for frame in stream:
    fft = frame['fft'][:plotSize]
    
    if 'mag' in plots:
        spec =  20.0 / scipy.log(10.0) * scipy.log(abs(fft) + 1e-7)

        pylab.subplot(len(plots), 1, plots.index('mag') + 1)
        pylab.gca().clear()
        
        pylab.gca().set_autoscale_on(False)
        pylab.gca().set_xlim([0, plotSize])
        pylab.gca().set_ylim([-100, 40])
        
        pylab.plot(spec)


    if 'phase' in plots:
        phase =  scipy.angle(fft)
        
        pylab.subplot(len(plots), 1, plots.index('phase') + 1)
        pylab.gca().clear()
        
        pylab.gca().set_autoscale_on(False)
        pylab.gca().set_xlim([0, plotSize])
        pylab.gca().set_ylim([-scipy.pi, scipy.pi])
        
        pylab.plot(phase)


    if 'peaks' in plots:
        fft = scipy.reshape(fft, (1, plotSize))
        peakLocs, peakMags, peakAng =  peaker.process( fft )
        
        pylab.subplot(len(plots), 1, plots.index('peaks') + 1)
        pylab.gca().clear()
        
        pylab.gca().set_autoscale_on(False)
        pylab.gca().set_xlim([0, plotSize])
        pylab.gca().set_ylim([-100, 40])
        
        pylab.stem(peakLocs, peakMags)


pylab.ioff()
