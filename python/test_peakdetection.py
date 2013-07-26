#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Create input
import scipy
import loudia

fundamental = 440.0
harmonics = 5

frameSize = 128
fftSize = 512
sampleRate = 8000

a_zeros = scipy.zeros( frameSize )
a_ones = scipy.ones( frameSize )
a_random = scipy.random.random( frameSize )
a_sine = scipy.cos(2 * scipy.pi * 440 * scipy.arange( frameSize ) / sampleRate + scipy.pi/4.0)

a_sine += (a_random - 0.5) * 1.0

# Loudia's solution # --------------------------------- #
window = loudia.Window(frameSize, loudia.Window.HAMMING)
fft = loudia.FFT(fftSize)
peaks = loudia.PeakDetection(5, loudia.PeakDetection.BYMAGNITUDE, 4)

r_sine_windowed = window.process(a_sine)
r_sine_fft = abs(fft.process(r_sine_windowed))

peaks.setCandidateCount( -1 )
peaks.setSortMethod( loudia.PeakDetectionComplex.BYMAGNITUDE )
r_sine_peakstart, r_sine_peakpos, r_sine_peakend, r_sine_peakmag = peaks.process(r_sine_fft)

peaks.setCandidateCount( -1 )
peaks.setSortMethod( loudia.PeakDetectionComplex.BYPOSITION )
r_sine_peakstart2, r_sine_peakpos2, r_sine_peakend2, r_sine_peakmag2 = peaks.process(r_sine_fft)

peaks.setCandidateCount( 10 )
peaks.setSortMethod( loudia.PeakDetectionComplex.BYMAGNITUDE )
r_sine_peakstart3, r_sine_peakpos3, r_sine_peakend3, r_sine_peakmag3 = peaks.process(r_sine_fft)

peaks.setCandidateCount( 10 )
peaks.setSortMethod( loudia.PeakDetectionComplex.BYPOSITION )
r_sine_peakstart4, r_sine_peakpos4, r_sine_peakend4, r_sine_peakmag4 = peaks.process(r_sine_fft)

# -------------------------------------------------------- #

import pylab

pylab.hold(True)
pylab.plot(r_sine_fft[0,:])

pylab.hold(True)
pylab.scatter(r_sine_peakpos[0,:], r_sine_peakmag[0,:], c = 'r', label = '-1, mag')
pylab.scatter(r_sine_peakstart[0,:], r_sine_fft[0, scipy.asarray(r_sine_peakstart[0,:], dtype='i')]*1.2, c = 'b', label = '-1, mag')
pylab.scatter(r_sine_peakend[0,:], r_sine_fft[0, scipy.asarray(r_sine_peakend[0,:], dtype='i')]*0.8, c = 'g', label = '-1, mag')

pylab.legend()


pylab.show()
