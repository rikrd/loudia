#!/usr/bin/env python

# Create input
import scipy
import loudia

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

# Loudia's solution # --------------------------------- #
window = loudia.Window(frameSize, loudia.Window.HAMMING)
fft = loudia.FFT(fftSize)
peaks = loudia.PeakDetectionComplex(5, 4)
peaksinterp = loudia.PeakInterpolationComplex()
trajs = loudia.PeakTracking(5, 4, 3)


r_sine_windowed = window.process(a_sine)
r_sine_mag = fft.process(r_sine_windowed)
r_sine_peakpos, r_sine_peakmag, r_sine_peakphase = peaks.process(r_sine_mag)
r_sine_peakipos, r_sine_peakimag, r_sine_peakphasei = peaksinterp.process(r_sine_mag, r_sine_peakpos, r_sine_peakmag, r_sine_peakphase)
r_sine_trajpos, r_sine_trajmag = trajs.process(r_sine_mag, r_sine_peakipos, r_sine_peakimag)
# -------------------------------------------------------- #

print r_sine_mag
print r_sine_peakpos
print r_sine_peakmag
print r_sine_trajpos, r_sine_trajmag

import pylab
pylab.hold(True)
pylab.plot(abs(r_sine_mag[0,:]))
pylab.hold(True)
pylab.scatter(r_sine_peakpos[0,:], r_sine_peakmag[0,:])

pylab.show()
