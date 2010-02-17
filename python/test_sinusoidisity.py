#!/usr/bin/env python

# Create input
import scipy
import loudia
import pylab
from scikits import audiolab

# Synthesize
filename = 'test-sinusoidisity.wav'

sampleRate = 44100.0    # samples / secs

duration = 10.0        # secs

frequency = 640.0      # Hz

vibratoAmplitude = frequency * (pow(2, 1.0/12.0) - 1.0) # Hz
vibratoFrequency = 5    # Hz

sampleCount = duration * sampleRate

noise = scipy.random.random(sampleCount)

t = scipy.arange(sampleCount) / sampleRate
phi_t = 2.0 * scipy.pi * frequency * t - vibratoAmplitude / vibratoFrequency * scipy.cos(2.0*scipy.pi*vibratoFrequency*t)
sine = scipy.sin(phi_t)

signal = sine + noise

# Writer sinus
audiolab.wavwrite(signal, filename, fs=sampleRate, enc='pcm16')

# Cut into frames
frameSize = 1024        # samples
hopSize = 512           # samples
fftSize = 2048          # samples

miniHopSize = 10        # samples


window = loudia.Window(frameSize, loudia.Window.BLACKMANHARRIS)
fft = loudia.FFT(fftSize, True)

def process(frame):
    spec = fft.process(window.process(frame))
    return scipy.absolute(spec), scipy.angle(spec)


cursor = 2*miniHopSize
while cursor < (sampleCount - frameSize):
    frameCurrent = signal[cursor:cursor+frameSize]
    framePast = signal[cursor-miniHopSize:cursor-miniHopSize+frameSize]
    framePast2 = signal[cursor-miniHopSize-miniHopSize:cursor-miniHopSize-miniHopSize+frameSize]

    cAbs, cAng = process(frameCurrent)
    pAbs, pAng = process(framePast)
    p2Abs, p2Ang = process(framePast2)

    cursor += hopSize
