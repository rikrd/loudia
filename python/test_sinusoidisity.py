#!/usr/bin/env python

# Create input
import scipy
import loudia
import pylab
from scikits import audiolab

filename = 'test-sinusoidisity.wav'

plot = False
frameSize = 1024        # samples
fftSize = 512           # samples
sampleRate = 44100.0    # samples / secs

duration = 100.0        # secs

frequency = 440.0      # Hz

vibratoAmplitude = frequency * (pow(2, 1.0/12.0) - 1.0) # Hz
vibratoFrequency = 5    # Hz

print vibratoAmplitude

sampleCount = duration * sampleRate

noise = scipy.random.random(sampleCount)

t = scipy.arange(sampleCount) / sampleRate
phi_t = 2.0 * scipy.pi * frequency * t - vibratoAmplitude / vibratoFrequency * scipy.cos(2.0*scipy.pi*vibratoFrequency*t)
sine = 0.5 * scipy.sin(phi_t)

# Writer sinus
audiolab.wavwrite(sine, filename, fs=sampleRate, enc='pcm16')
