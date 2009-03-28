#!/usr/bin/env python

# Create input
import scipy
import ricaudio

fundamental = 1110.0
harmonics = 10

frameSize = 128
fftSize = 4096
samplerate = 44100.0

a_zeros = scipy.zeros( frameSize )
a_ones = scipy.ones( frameSize )
a_random = scipy.random.random( frameSize )
a_sine = scipy.cos(2 * scipy.pi * 440 * scipy.arange( frameSize ) / samplerate + scipy.pi/4.0)

a_sine += (a_random - 0.5) * 1.0

# Ricaudio's solution # --------------------------------- #
m = ricaudio.SpectralReassignment(frameSize, fftSize, samplerate, ricaudio.Window.BLACKMANHARRIS)

r_zeros_fft, r_zeros_time, r_zeros_freq = m.process(a_zeros)
r_ones_fft, r_ones_time, r_ones_freq = m.process(a_ones)
r_random_fft, r_random_time, r_random_freq = m.process(a_random)
r_sine_fft, r_sine_time, r_sine_freq = m.process(a_sine)
# -------------------------------------------------------- #


# Scipy solution # ---------------------------------- #
s_zeros_fft = scipy.fft(a_zeros, fftSize)
s_ones_fft = scipy.fft(a_ones, fftSize)
s_random_fft = scipy.fft(a_random, fftSize)
s_sine_fft = scipy.fft(a_sine, fftSize)
# -------------------------------------------------------- #

print "FFT:"
print r_zeros_fft
print r_ones_fft
print r_random_fft
print r_sine_fft

print "Time:"
print r_zeros_time
print r_ones_time
print r_random_time
print r_sine_time

print "Freq:"
print r_zeros_freq
print r_ones_freq
print r_random_freq
print r_sine_freq

print scipy.allclose(r_zeros_fft, s_zeros_fft)
print scipy.allclose(r_ones_fft, s_ones_fft)
print scipy.allclose(r_random_fft, s_random_fft)
print scipy.allclose(r_sine_fft, s_sine_fft)

r_abs = scipy.absolute(r_sine_fft).T
r_ang = scipy.angle(r_sine_fft).T
r_max = max(r_abs)

s_abs = scipy.absolute(s_sine_fft).T
s_ang = scipy.angle(s_sine_fft).T
s_max = max(s_abs)

import pylab
pylab.subplot(311)
pylab.plot(a_sine.T)

pylab.subplot(312)
pylab.hold(True)
pylab.plot(r_abs, label = 'Ricaudio')
pylab.plot(s_abs, label = 'Scipy')
pylab.legend()

pylab.subplot(313)
pylab.hold(True)
pylab.plot(r_abs*r_ang/r_max, label = 'Ricaudio')
pylab.plot(s_abs*s_ang/s_max, label = 'Scipy')
pylab.legend()

pylab.figure()
pylab.subplot(211)
pylab.plot(r_sine_time.T, label = 'Time Reqssignment')
pylab.legend()

pylab.subplot(212)
pylab.plot(r_sine_freq.T, label = 'Frequency Reqssignment')
pylab.legend()
pylab.show()
