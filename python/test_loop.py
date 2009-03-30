#!/usr/bin/env python

# Create input
import scipy
import loudia

frameSize = 121
fftSize = 512
samplerate = 8000


# Loudia's solution # --------------------------------- #
w = loudia.Window(frameSize, 0)
m = loudia.FFT(fftSize, True)

for i in range(100000):
    a_random = scipy.random.random(frameSize)
    r_random = m.process(w.process(a_random))
# -------------------------------------------------------- #
