#!/usr/bin/env python

# Create input
import scipy
import ricaudio

fundamental = 440.0
harmonics = 10

frameSize = 1024
frameStep = 1
beginSize = frameSize
signalSize = 2*frameSize
fftSize = frameSize * (2**2)
samplerate = 44100

sine = scipy.zeros((signalSize, 1), dtype = 'f4')
for i in range(harmonics):
    sine[beginSize:, 0] += scipy.array(scipy.cos(2 * scipy.pi * i * fundamental * scipy.arange(signalSize - beginSize) / samplerate), dtype='f4')

sinenoise = scipy.array(sine + (scipy.random.random((signalSize, 1)) - 0.5) * 0.00125, dtype='f4')

# Ricaudio's solution # --------------------------------- #
bandwidth = 8
peakCount = 10
window = ricaudio.Window(frameSize, ricaudio.Window.BLACKMANHARRIS)
fft = ricaudio.FFT(frameSize, fftSize)
peaks = ricaudio.PeakDetect(peakCount, bandwidth)
cogs = ricaudio.PeakCOG(fftSize, bandwidth)

peakCogs = []
for i in range((signalSize - frameSize)/frameStep):
    a_sine = sine[i*frameStep:(i*frameStep+frameSize)].T
  
    r_sine_windowed = window.process(a_sine)
    r_sine_fft = fft.process(r_sine_windowed)
    r_sine_fft = r_sine_fft[:,:scipy.ceil(r_sine_fft.shape[1]/2)]
    r_sine_peakpos, r_sine_peakmag, r_sine_peakphase = peaks.process(r_sine_fft)
    r_sine_peakcogs = cogs.process(r_sine_fft, r_sine_peakpos)

    peakCogs.append(r_sine_peakcogs[0,:5])
# -------------------------------------------------------- #

peakCogs = scipy.array(peakCogs)
print peakCogs

import pylab
pylab.figure()
pylab.plot(peakCogs)
pylab.show()


