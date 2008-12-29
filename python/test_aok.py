#!/usr/bin/env python

# Create input
import scipy
windowSize = 300.0
hopSize = 16000.0
fftSize = 512
normVolume = 3.0

a1 = scipy.array(scipy.zeros((1, 512)), dtype='f4')
a2 = scipy.array(scipy.ones((1, 512)), dtype='f4')

# CRicaudio's solution # --------------------------------- #
import ricaudio
m = ricaudio.AOK(windowSize, hopSize, fftSize, normVolume)

b1 = m.process(a1)
b2 = m.process(a2)
# -------------------------------------------------------- #

print b1
print b2
