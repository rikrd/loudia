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
peaks = loudia.PeakDetectionComplex(5, loudia.PeakDetectionComplex.BYMAGNITUDE, 4)

r_sine_windowed = window.process(a_sine)
r_sine_fft = fft.process(r_sine_windowed)

peaks.setCandidateCount( -1 )
peaks.setSortMethod( loudia.PeakDetectionComplex.BYMAGNITUDE )
r_sine_peakpos, r_sine_peakmag, r_sine_peakphase = peaks.process(r_sine_fft)

peaks.setCandidateCount( -1 )
peaks.setSortMethod( loudia.PeakDetectionComplex.BYPOSITION )
r_sine_peakpos2, r_sine_peakmag2, r_sine_peakphase2 = peaks.process(r_sine_fft)

peaks.setCandidateCount( 10 )
peaks.setSortMethod( loudia.PeakDetectionComplex.BYMAGNITUDE )
r_sine_peakpos3, r_sine_peakmag3, r_sine_peakphase3 = peaks.process(r_sine_fft)

peaks.setCandidateCount( 10 )
peaks.setSortMethod( loudia.PeakDetectionComplex.BYPOSITION )
r_sine_peakpos4, r_sine_peakmag4, r_sine_peakphase4 = peaks.process(r_sine_fft)

# -------------------------------------------------------- #

import pylab
pylab.subplot(211)
pylab.hold(True)
pylab.plot(abs(r_sine_fft[0,:]))
pylab.hold(True)
pylab.scatter(r_sine_peakpos[0,:], r_sine_peakmag[0,:], c = 'r', label = '-1, mag')
pylab.scatter(r_sine_peakpos2[0,:], r_sine_peakmag2[0,:], c = 'b', label = '-1, pos')
pylab.scatter(r_sine_peakpos3[0,:], r_sine_peakmag3[0,:], c = 'k', label = '10, mag')
pylab.scatter(r_sine_peakpos4[0,:], r_sine_peakmag4[0,:], c = 'g', label = '10, pos')

pylab.legend()

pylab.subplot(212)
pylab.hold(True)
pylab.plot(scipy.angle(r_sine_fft[0,:]))
pylab.hold(True)
pylab.scatter(r_sine_peakpos[0,:], r_sine_peakphase[0,:], c = 'r', label = '-1, mag')
pylab.scatter(r_sine_peakpos2[0,:], r_sine_peakphase2[0,:], c = 'b', label = '-1, pos')
pylab.scatter(r_sine_peakpos3[0,:], r_sine_peakphase3[0,:], c = 'k', label = '10, mag')
pylab.scatter(r_sine_peakpos4[0,:], r_sine_peakphase4[0,:], c = 'g', label = '10, pos')

pylab.legend()

pylab.show()
