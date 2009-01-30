#!/usr/bin/env python

import scipy
import pylab
import ricaudio
from sepel.inputs import pyricaudio
import sys

filename = sys.argv[1]

frameSize = 256 
frameStep = 256

frameSizeTime = frameSize / 44100.0
frameStepTime = frameStep / 44100.0

fftSize = frameSize * (2**3)

analysisLimit = 1000.0

plotSize = fftSize / 4

# Creation of the pipeline        
stream = pyricaudio.sndfilereader({'filename': filename,
                                   'windowSizeInTime': frameSizeTime,
                                   'windowStepInTime': frameStepTime,
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

pylab.figure()
for ind, frame in enumerate(stream):
    if ind == 34:
        spectrum =  20.0 / scipy.log( 10.0 ) * scipy.log( abs( frame['fft'][:plotSize] ) + 1e-7)
        pylab.plot(spectrum)
        pylab.show()
