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
fft = ricaudio.FFT(fftSize)
peaks = ricaudio.PeakDetectComplex(fftSize / 3)
peaksi = ricaudio.PeakInterpolateComplex()

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
