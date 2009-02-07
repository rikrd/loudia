#!/usr/bin/env python

from common import *

plot = True

filename = '/home/rmarxer/dev/data/onsets/pitchedphrases/Strings/Picked-Plucked-Ham/Piano/piano1.wav'
windowSize = 1024
windowHop = 512
windowType = ricaudio.Window.HAMMING
fftSize = 2048
halfSize = 1025
zeroPhase = True

components = 10
pastFrames = 10
pastCoeff = 0.2

w = ricaudio.Window(windowSize, windowType)
f = ricaudio.FFT(windowSize, fftSize, zeroPhase)
d = ricaudio.INMF(halfSize, components, pastFrames, pastCoeff, 15, 1, 1e-17)

framer, sr, nframes, nchannels, loader = get_framer_audio(filename, windowSize, windowHop)

components = []
gains = []
a = []
for frame in framer:
    windowed = w.process(frame)
    fft = f.process(windowed)
    fft = abs(fft)
    c, g = d.process(fft)

    a.append(20.0 * scipy.log10(fft[0,:]+1e-9))
    gains.append(g)
    components.append(20.0 * scipy.log10(c[0,:]+1e-9))

gains = overlapadder(gains, 1, 1)
components = components[-1]
a = scipy.array(a)

if plot:
    pylab.figure()
    pylab.imshow(scipy.flipud(a.T))
    pylab.title("Input")
    
    pylab.figure()
    pylab.plot(components.T)
    pylab.title("Components")
    
    pylab.figure()
    pylab.plot(gains)
    pylab.title("Gains")
    
    pylab.show()
