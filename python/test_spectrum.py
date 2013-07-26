#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

plotSize = fftSize

stream, sampleRate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)

windower = loudia.Window( frameSize, loudia.Window.HAMMING )
ffter = loudia.FFT( fftSize )

pylab.figure()
for ind, frame in enumerate(stream):
    if ind == 34:
        fft = ffter.process( windower.process( frame ) )
        peaks = detectPeaks( abs(fft)[0], tol = 0)
        starts, pos, mags = zip(*peaks)
        #spec =  loudia.magToDb( abs( fft ) )
        spec =  abs( fft ) 
        starts = scipy.array(starts)
        pos = scipy.array(pos)
        
        pylab.plot( spec[0, :plotSize] )
        pylab.stem( starts[starts < plotSize], [spec[0, s] for s in starts[ starts < plotSize ]], linefmt=' ', markerfmt='go' )
        pylab.stem( pos[pos < plotSize], [spec[0, s] for s in pos[ pos < plotSize ]] )
        pylab.show()
