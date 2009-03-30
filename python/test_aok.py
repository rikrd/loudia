#!/usr/bin/env python

# Create input
import scipy
windowSize = 300.0
hopSize = 16000.0
fftSize = 512
normVolume = 3.0

a1 = scipy.zeros(512)
a2 = scipy.ones(512)

# Loudia's solution # --------------------------------- #
import loudia
m = loudia.AOK(windowSize, hopSize, fftSize, normVolume)

b1 = m.process(a1)
b2 = m.process(a2)
# -------------------------------------------------------- #

print b1
print b2
