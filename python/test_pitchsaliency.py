#!/usr/bin/env python

import ricaudio
from sepel.inputs import pyricaudio
import pylab
import os, sys, wave
import scipy
from common import *

interactivePlot = False
plot = True

filename = sys.argv[1]

# Samplerate of the file
wavfile = wave.open(filename,'r')
samplerate = float(wavfile.getframerate())
wavfile.close()

frameSize = 2048 
frameStep = 1024

frameSizeTime = frameSize / 44100.0
frameStepTime = frameStep / 44100.0

fftSize = 4096
plotSize = fftSize / 4

bandwidth = 4 * fftSize/frameSize
analysisLimit = scipy.inf

# Creation of the pipeline        
stream = pyricaudio.sndfilereader({'filename': filename,
                                   'windowSizeInTime': frameSizeTime,
                                   'windowStepInTime': frameStepTime,
                                   'encodingKey': 'encoding',
                                   'channelCountKey': 'channelCount',
                                   'samplesOriginalKey': 'samples',
                                   'samplesKey': 'samplesMono',
                                   'samplerateKey': 'samplerate',
                                   'timestampBeginKey': 'timestampBegin',
                                   'timestampEndKey': 'timestampEnd',
                                   'limit':analysisLimit})

stream = pyricaudio.window_ricaudio(stream, {'inputKey': 'samplesMono',
                                             'outputKey': 'windowed',
                                             'windowType': 'hamming'})

stream = pyricaudio.fft_ricaudio(stream, {'inputKey': 'windowed',
                                          'outputKey': 'fft',
                                          'zeroPhase': True,
                                          'fftLength': fftSize})

freqPrec = 0.001
deltaPeriod = 2
numHarmonics = 10

whitening = ricaudio.SpectralWhitening(fftSize, 50.0, 6000.0, samplerate)
pitchSaliency = ricaudio.PitchSaliency(fftSize, 50.0, 2100.0, samplerate, freqPrec, numHarmonics)

specs = []
wspecs = []
pitches = []
saliencies = []

if interactivePlot:
    pylab.ion()
    pylab.figure()
    pylab.title('Interactive plot of the FFT vs LPC frequency response')
    pylab.gca().set_ylim([-100, 40])
    pylab.gca().set_autoscale_on(False)

for frame in stream:
    spec = scipy.array(abs(frame['fft']), dtype = scipy.float32)

    wspec = whitening.process( spec )
    pitch, saliency = pitchSaliency.process( wspec )
    
    if interactivePlot:
        pylab.subplot(211)
        pylab.hold(False)
        pylab.plot( spec )
        
        pylab.subplot(212)
        pylab.hold(False)
        pylab.plot(result[0,:], label = 'Noise Suppressed Spectrum')
        
    specs.append( spec )
    wspecs.append( wspec )
    pitches.append( pitch )
    saliencies.append( saliency )

if interactivePlot:
    pylab.ioff()


specs = scipy.array( specs )
wspecs = scipy.array( wspecs )
pitches = scipy.array( pitches )[:, 0]
saliencies = scipy.array( saliencies )[:, 0]
frameCount = specs.shape[0] - 1

if plot:
    pylab.figure()
    pylab.subplot(211)
    pylab.plot( pitches[:,0] )

    pylab.subplot(212)
    pylab.plot( saliencies[:,0] )
    pylab.show()

print pitches
