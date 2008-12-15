#!/usr/bin/env python

# Create input
import scipy
import ricaudio

fundamental = 440.0
harmonics = 5

frameSize = 256
fftSize = 512
samplerate = 8000

a_zeros = scipy.array(scipy.zeros((1, frameSize)), dtype='f4')
a_ones = scipy.array(scipy.ones((1, frameSize)), dtype='f4')
a_random = scipy.array(scipy.random.random((1, frameSize)), dtype='f4')

a_sine = a_zeros
for i in range(harmonics):
    a_sine[0, :] += scipy.array(scipy.cos(2 * scipy.pi * i * fundamental * scipy.arange(frameSize) / samplerate), dtype='f4')

a_sine += (a_random - 0.5) * 1.0

# Ricaudio's solution # --------------------------------- #
m = ricaudio.SpectralReassignment(frameSize, fftSize, samplerate, 2)

r_zeros = m.process(a_zeros)
r_ones = m.process(a_ones)
r_random = m.process(a_random)
r_sine = m.process(a_sine)
# -------------------------------------------------------- #


# Scipy solution # ---------------------------------- #
s_zeros = scipy.fft(a_zeros, fftSize)
s_ones = scipy.fft(a_ones, fftSize)
s_random = scipy.fft(a_random, fftSize)
s_sine = scipy.fft(a_sine, fftSize)
# -------------------------------------------------------- #


print scipy.allclose(r_zeros, s_zeros)
print scipy.allclose(r_ones, s_ones)
print scipy.allclose(r_random, s_random)
print scipy.allclose(r_sine, s_sine)

r_abs = scipy.absolute(r_sine).T
r_ang = scipy.angle(r_sine).T
r_max = max(r_abs)

s_abs = scipy.absolute(s_sine).T
s_ang = scipy.angle(s_sine).T
s_max = max(s_abs)

import pylab
pylab.subplot(311)
pylab.plot(a_sine.T)

pylab.subplot(312)
pylab.hold(True)
pylab.plot(r_abs, label = 'Ricaudio')
pylab.plot(s_abs, label = 'Scipy')

pylab.subplot(313)
pylab.hold(True)
pylab.plot(r_abs*r_ang/r_max, label = 'Ricaudio')
pylab.plot(s_abs*s_ang/s_max, label = 'Scipy')

pylab.legend()
pylab.show()
