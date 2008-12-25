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
window = ricaudio.Window(frameSize, ricaudio.Window.HAMMING)
fft = ricaudio.FFT(frameSize, fftSize)
peaks = ricaudio.SpectralPeaks(-1)

r_sine_windowed = window.process(a_sine)
r_sine_mag = abs(fft.process(r_sine_windowed))
r_sine_pos, r_sine_mag = peaks.process(r_sine_mag)
# -------------------------------------------------------- #

print r_sine_pos
print r_sine_mag

