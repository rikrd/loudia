#!/usr/bin/env python

import loudia
from common import *
import pylab
import os, sys, wave
import scipy

interactivePlot = False
plot = True

filename = sys.argv[1]

frameSize = 4096
frameStep = 1024

fftSize = 4096
plotSize = fftSize / 4

bandwidth = 4 * fftSize/frameSize

minPeakWidth = 8

peakCandidateCount = 4

stream, samplerate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)

ffter = loudia.FFT( fftSize )
windower = loudia.Window( frameSize, loudia.Window.BLACKMANHARRIS )
whitening = loudia.SpectralWhitening(fftSize, 50.0, 6000.0, samplerate)
pitchACF = loudia.PitchACF(fftSize, samplerate, minPeakWidth, peakCandidateCount)
acorr = loudia.Autocorrelation(fftSize/2+1, fftSize/2+1)

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
    fft = ffter.process( windower.process( frame ) )
    spec =  loudia.magToDb( abs( fft ) )
    
    wspec = whitening.process( spec )
    pitch, saliency = pitchACF.process( wspec )
    
    if interactivePlot:
        pylab.subplot(211)
        pylab.hold(False)
        pylab.plot( wspec[0, :plotSize] )

        acorred = acorr.process( wspec )
        
        pylab.subplot(212)
        pylab.hold(False)
        pylab.plot(acorred[0,:plotSize], label = 'Noise Suppressed Spectrum')
        pylab.hold(True)
        pylab.stem( pitch/samplerate*fftSize, saliency )
        
    specs.append( spec[0, :] )
    wspecs.append( wspec[0, :] )
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
    saliencies[ pitches < 0 ] = scipy.nan
    pitches[ pitches < 0 ] = scipy.nan
    
    # Get the onsets
    annotation = os.path.splitext(filename)[0] + '.onset_annotated'
    onsets = get_onsets(annotation, frameStep, samplerate)
    
    pylab.figure()
    pylab.subplot(211)
    pylab.plot( pitches[:,0] )

    draw_onsets(onsets)

    pylab.subplot(212)
    pylab.plot( saliencies[:,0] )

    draw_onsets(onsets)
    
    pylab.show()

print pitches
