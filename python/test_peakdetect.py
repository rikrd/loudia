#!/usr/bin/env python

# Create input
import scipy
import ricaudio

fundamental = 440.0
harmonics = 5

frameSize = 128
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
peaks = ricaudio.PeakDetect(5, 4)

r_sine_windowed = window.process(a_sine)
r_sine_fft = fft.process(r_sine_windowed)
r_sine_fft = r_sine_fft[:,:scipy.ceil(r_sine_fft.shape[1]/2)]
r_sine_peakpos, r_sine_peakmag, r_sine_peakphase = peaks.process(r_sine_fft)
# -------------------------------------------------------- #

print r_sine_fft
print r_sine_peakpos
print r_sine_peakmag
print r_sine_peakphase

import pylab
pylab.subplot(211)
pylab.hold(True)
pylab.plot(abs(r_sine_fft[0,:]))
pylab.hold(True)
pylab.scatter(r_sine_peakpos[0,:], r_sine_peakmag[0,:])

pylab.subplot(212)
pylab.hold(True)
pylab.plot(scipy.angle(r_sine_fft[0,:]))
pylab.hold(True)
pylab.scatter(r_sine_peakpos[0,:], r_sine_peakphase[0,:])


pylab.show()
