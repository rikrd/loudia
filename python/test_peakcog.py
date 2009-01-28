#!/usr/bin/env python

# Create input
import pylab
import scipy
import ricaudio

frameSize = 1024
frameStep = 64
beginSize = frameSize
signalSize = 2*frameSize
fftSize = frameSize * (2**2)
samplerate = 44100

fundamental = 0.100 * samplerate
harmonics = 20

plotInteractive = False

sine = scipy.zeros((signalSize, 1), dtype = 'f4')
for i in range(harmonics):
    sine[beginSize:, 0] += scipy.array(scipy.cos(2 * scipy.pi * (i+1) * fundamental * scipy.arange(signalSize - beginSize) / samplerate), dtype='f4')

sinenoise = scipy.array(sine + (scipy.random.random((signalSize, 1)) - 0.5), dtype='f4')

# Ricaudio's solution # --------------------------------- #
bandwidth = 8 * fftSize / frameSize
peakCount = 4
window = ricaudio.Window(frameSize, ricaudio.Window.BLACKMANHARRIS)
fft = ricaudio.FFT(frameSize, fftSize)
peaks = ricaudio.PeakDetect(peakCount, bandwidth)
cogs = ricaudio.PeakCOG(fftSize, bandwidth)

peakCogs = []
peakPos = []
if plotInteractive:
    pylab.ion()
    pylab.hold(False)
    
for i in range((signalSize - frameSize)/frameStep):
    a_sine = sinenoise[i*frameStep:(i*frameStep+frameSize)].T
  
    r_sine_windowed = window.process(a_sine)
    r_sine_fft = fft.process(r_sine_windowed)
    r_sine_fft = r_sine_fft[:,:scipy.ceil(r_sine_fft.shape[1]/2)]
    r_sine_peakpos, r_sine_peakmag, r_sine_peakphase = peaks.process(r_sine_fft)
    r_sine_peakcogs = cogs.process(r_sine_fft, r_sine_peakpos)

    if plotInteractive:
        pylab.clf()
        pylab.hold(True)
        pylab.plot(abs(r_sine_fft[0,:]))
        pylab.stem(r_sine_peakpos[0,:], abs(r_sine_fft[0,:])[scipy.array(r_sine_peakpos[0,:], dtype = 'i4')])
        
        pylab.hold(False)
        
    peakPos.append(r_sine_peakpos[0,:])
    peakCogs.append(r_sine_peakcogs[0,:])
# -------------------------------------------------------- #
if plotInteractive:
    pylab.ioff()
    
peakCogs = scipy.array(peakCogs)
peakPos = scipy.array(peakPos)
print peakPos[:10,:]

pylab.figure()
pylab.plot(peakCogs)
pylab.show()


