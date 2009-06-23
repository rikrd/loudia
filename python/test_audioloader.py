#!/usr/bin/env python

# Create input
import scipy
import loudia
from common import *
import sys
import pylab

interactivePlotting = False

atol = 1e-5
rtol = 1e-5

filename = sys.argv[1]

frameSize = 1024 
frameStep = 1024

stream, sampleRate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)

# Loudia's solution # --------------------------------- #
loader = loudia.AudioLoader()
loader.setFilename(filename, False)
loader.setFrameSize(frameSize, False)
loader.setChannel(loudia.AudioLoader.MIX, False)
loader.setup()

if interactivePlotting:
    pylab.ion()
    
for audiolabFrame in stream[:10]:
    loudiaFrame = loader.process()[:, 0]

    if interactivePlotting:
        pylab.subplot(211)
        pylab.gca().clear()
        pylab.plot(loudiaFrame)
        pylab.subplot(212)
        pylab.gca().clear()
        pylab.plot(audiolabFrame)
    
    print scipy.allclose(loudiaFrame, audiolabFrame, atol = atol, rtol = rtol)

if interactivePlotting:
    pylab.ioff()
