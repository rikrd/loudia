#!/usr/bin/env python

# Create input
import scipy
import ricaudio

frameSize = 121
fftSize = 512
samplerate = 8000


# Ricaudio's solution # --------------------------------- #
w = ricaudio.Window(frameSize, 0)
m = ricaudio.FFT(fftSize, True)

for i in range(100000):
    a_random = scipy.random.random(frameSize)
    r_random = m.process(w.process(a_random))
# -------------------------------------------------------- #
