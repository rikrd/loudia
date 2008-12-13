#!/usr/bin/env python

# Create input
import scipy
import ricaudio

frameSize = 64
fftSize = 128
samplerate = 8000

a_zeros = scipy.array(scipy.zeros((1, frameSize)), dtype='f4')
a_ones = scipy.array(scipy.ones((1, frameSize)), dtype='f4')
a_random = scipy.array(scipy.random.random((1, frameSize)), dtype='f4')
a_sine = scipy.array(scipy.cos(2 * scipy.pi * 440 * scipy.arange(frameSize) / samplerate), dtype='f4')
a_sine = a_sine.reshape((1, a_sine.shape[0]))

# Ricaudio's solution # --------------------------------- #
m = ricaudio.FFT(frameSize, fftSize, True)

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



import pylab
pylab.hold(True)
pylab.subplot(211)
pylab.plot(scipy.absolute(r_sine).T, label = 'Ricaudio')
pylab.plot(scipy.absolute(s_sine), label = 'Scipy')

pylab.subplot(212)
pylab.plot(scipy.angle(r_sine).T, label = 'Ricaudio')
pylab.plot(scipy.angle(s_sine), label = 'Scipy')

pylab.show()
