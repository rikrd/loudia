#!/usr/bin/env python

# Create input
import scipy
import ricaudio

fundamental = 440.0
harmonics = 1

frameSize = 1024
frameStep = 128
beginSize = frameSize
signalSize = 2*frameSize
fftSize = frameSize * (2**2)
samplerate = 44100

sine = scipy.zeros((signalSize, 1), dtype = 'f4')
for i in range(harmonics):
    print scipy.arange(signalSize - beginSize).shape
    sine[beginSize:, 0] += scipy.array(scipy.cos(2 * scipy.pi * i * fundamental * scipy.arange(signalSize - beginSize) / samplerate), dtype='f4')

sinenoise = sine + (scipy.random.random((signalSize, 1)) - 0.5) * 1.0

# Ricaudio's solution # --------------------------------- #
bandwidth = 8
window = ricaudio.Window(frameSize, ricaudio.Window.BLACKMANHARRIS)
fft = ricaudio.FFT(frameSize, fftSize)
peaks = ricaudio.PeakDetect(fftSize / 3, bandwidth)
cogs = ricaudio.PeakCOG(fftSize, bandwidth)

peakCogs = []
for i in range((signalSize - frameSize)/frameStep):
    a_sine = sine[i*frameStep:(i*frameStep+frameSize)].T
  
    r_sine_windowed = window.process(a_sine)
    r_sine_fft = fft.process(r_sine_windowed)
    r_sine_fft = r_sine_fft[:,:scipy.ceil(r_sine_fft.shape[1]/2)]
    r_sine_peakpos, r_sine_peakmag, r_sine_peakphase = peaks.process(r_sine_fft)
    r_sine_peakcogs = cogs.process(r_sine_fft, r_sine_peakpos)

    peakCogs.append(r_sine_peakcogs[0,:])
# -------------------------------------------------------- #

peakCogs = scipy.array(peakCogs)
import pylab
pylab.figure()
pylab.plot(peakCogs)
pylab.show()


