#!/usr/bin/env python

# Create input
import scipy
import ricaudio

frameSize = 256
fftSize = 512

a_zeros = scipy.array(scipy.zeros((1, frameSize)), dtype='f4')
a_ones = scipy.array(scipy.ones((1, frameSize)), dtype='f4')
a_random = scipy.array(scipy.random.random((1, frameSize)), dtype='f4')

# Ricaudio's solution # --------------------------------- #
m = ricaudio.FFT(frameSize, fftSize)

r_zeros = m.process(a_zeros)
r_ones = m.process(a_ones)
r_random = m.process(a_random)
# -------------------------------------------------------- #


# Scipy solution # ---------------------------------- #
s_zeros = scipy.fft(a_zeros, fftSize)
s_ones = scipy.fft(a_ones, fftSize)
s_random = scipy.fft(a_random, fftSize)
# -------------------------------------------------------- #

print scipy.allclose(r_zeros, s_zeros)
print scipy.allclose(r_ones, s_ones)
print scipy.allclose(r_random, s_random)
