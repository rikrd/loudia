#!/usr/bin/env python

import scipy
from common import *
import pylab
import loudia
import sys

filename = sys.argv[1]

frameSize = 256 
frameStep = 256

fftSize = frameSize * (2**3)

analysisLimit = 1000.0

plotSize = fftSize / 4

stream, samplerate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)

windower = loudia.Window( frameSize, loudia.Window.HAMMING )
ffter = loudia.FFT( fftSize )

pylab.figure()
for ind, frame in enumerate(stream):
    if ind == 34:
        fft = ffter.process( windower.process( frame ) )
        spec =  loudia.magToDb( abs( fft ) )

        pylab.plot( spec[0, :plotSize] )
        pylab.show()
