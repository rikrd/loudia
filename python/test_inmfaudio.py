#!/usr/bin/env python

from common import *
import sys

plot = True

filename = '/home/rmarxer/dev/data/onsets/pitchedphrases/Strings/Picked-Plucked-Ham/Piano/piano1.wav'

if len(sys.argv) >= 2:
    filename = sys.argv[1]

windowSize = 1024
windowHop = 256
windowType = ricaudio.Window.BLACKMANHARRIS
fftSize = 2048
halfSize = (fftSize / 2) + 1
zeroPhase = True

components = 2
pastFrames = 10
pastCoeff = 0.2

w = ricaudio.Window(windowSize, windowType)
f = ricaudio.FFT(windowSize, fftSize, zeroPhase)
d = ricaudio.INMF(halfSize, components, pastFrames, pastCoeff, 15, 1, 1e-17)

framer, sr, nframes, nchannels, loader = get_framer_audio(filename, windowSize, windowHop)

onsets = get_onsets(filename, windowHop, sr)

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
nwindows = a.shape[0]
if plot:
    pylab.figure()
    pylab.subplot(211)
    pylab.imshow(scipy.flipud(a.T), aspect = 'auto')
    draw_onsets(onsets)
    pylab.title("Input")
    pylab.gca().set_xlim([0, nwindows - 1])


    pylab.subplot(212)
    pylab.plot(gains)
    draw_onsets(onsets)
    pylab.title("Gains")
    pylab.gca().set_xlim([0, nwindows - 1])


    pylab.figure()
    pylab.plot(components.T)
    pylab.title("Components")

    pylab.show()
