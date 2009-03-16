#!/usr/bin/env python

# Create input
import scipy
import ricaudio

fundamental = 440.0
harmonics = 5

frameSize = 128
fftSize = 512
samplerate = 8000

a_zeros = scipy.zeros( frameSize )
a_ones = scipy.ones( frameSize )
a_random = scipy.random.random( frameSize )
a_sine = scipy.cos(2 * scipy.pi * 440 * scipy.arange( frameSize ) / samplerate + scipy.pi/4.0)

a_sine += (a_random - 0.5) * 1.0

# Ricaudio's solution # --------------------------------- #
window = ricaudio.Window(frameSize, ricaudio.Window.HAMMING)
fft = ricaudio.FFT(fftSize)
peaks = ricaudio.PeakPick(fftSize / 3)

r_sine_windowed = window.process(a_sine)
r_sine_mag = fft.process(r_sine_windowed)
r_sine_mag = r_sine_mag[:,:scipy.ceil(r_sine_mag.shape[1]/2)]
r_sine_peakpos, r_sine_peakmag = peaks.process(r_sine_mag)
# -------------------------------------------------------- #

print r_sine_mag
print r_sine_peakpos
print r_sine_peakmag

import pylab
pylab.hold(True)
pylab.plot(abs(r_sine_mag[0,:]))
pylab.hold(True)
pylab.scatter(r_sine_peakpos[0,:], r_sine_peakmag[0,:])

pylab.show()
