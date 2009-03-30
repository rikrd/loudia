#!/usr/bin/env python

# Create input
import scipy
import loudia

fundamental = 440.0
harmonics = 5

frameSize = 256
fftSize = 512
samplerate = 8000

a_zeros = scipy.zeros( frameSize )
a_ones = scipy.ones( frameSize )
a_random = scipy.random.random( frameSize )
a_sine = scipy.cos(2 * scipy.pi * 440 * scipy.arange( frameSize ) / samplerate + scipy.pi/4.0)

a_sine += (a_random - 0.5) * 1.0

# Ricaudio's solution # --------------------------------- #
window = loudia.Window(frameSize, loudia.Window.HAMMING)
fft = loudia.FFT(fftSize)
peaks = loudia.PeakDetectionComplex(fftSize / 3)
peaksi = loudia.PeakInterpolationComplex()

r_sine_windowed = window.process(a_sine)
r_sine_fft = fft.process(r_sine_windowed)
r_sine_fft = r_sine_fft[:,:scipy.ceil(r_sine_fft.shape[1]/2)]
r_sine_peakpos, r_sine_peakmag, r_sine_peakphase = peaks.process(r_sine_fft)
r_sine_peakposi, r_sine_peakmagi, r_sine_peakphasei = peaksi.process(r_sine_fft, r_sine_peakpos, r_sine_peakmag, r_sine_peakphase)

#r_sine_peakphasei[r_sine_peakphasei != -1] = ((r_sine_peakphasei[r_sine_peakphasei != -1] + scipy.pi) % (-2*scipy.pi)) + (scipy.pi)
# -------------------------------------------------------- #

print r_sine_fft
print r_sine_peakpos
print r_sine_peakmag
print r_sine_peakposi
print r_sine_peakmagi
print r_sine_peakphase
print r_sine_peakphasei

import pylab
pylab.subplot(211)
pylab.hold(True)
pylab.plot(20.0 * scipy.log10(abs(r_sine_fft[0,:])))
pylab.hold(True)
pylab.scatter(r_sine_peakpos[0,:], 20.0 * scipy.log10(r_sine_peakmag[0,:]))
pylab.scatter(r_sine_peakposi[0,:], r_sine_peakmagi[0,:], c='r')

pylab.subplot(212)
pylab.hold(True)
pylab.plot(scipy.angle(r_sine_fft[0,:]))
pylab.hold(True)
pylab.scatter(r_sine_peakpos[0,:], r_sine_peakphase[0,:])
pylab.scatter(r_sine_peakposi[0,:], r_sine_peakphasei[0,:], c='r')

pylab.show()
