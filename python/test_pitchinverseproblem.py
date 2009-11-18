#!/usr/bin/env python

import loudia
from common import *
import pylab
import os, sys, wave
import scipy


interactivePlot = True
plot = True

filename = sys.argv[1]

frameSize = 4096 
frameStep = 2048

fftSize = 8192

saliencyThreshold = 20

plotSize = fftSize / 8

stream, sampleRate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)


peakBandwidth = 4
peakCandidateCount = 4
numMaxPitches = 1
numHarmonics = 80
numCandidates = 200

windower = loudia.Window( frameSize,
                          loudia.Window.BLACKMANHARRIS )

ffter = loudia.FFT( fftSize )

whitening = loudia.SpectralWhitening(fftSize,
                                     100.0, 4100.0,
                                     sampleRate)

pitchInverseProblem = loudia.PitchInverseProblem(fftSize,
                                                 50.0, 2000.0,
                                                 sampleRate,
                                                 numMaxPitches,
                                                 numHarmonics,
                                                 numCandidates,
                                                 peakBandwidth)

a = pitchInverseProblem.projectionMatrix()

specs = []
wspecs = []
pitches = []
saliencies = []
freqss = []

frameNumber = 200

if interactivePlot:
    pylab.ion()
    pylab.figure()
    pylab.title('Interactive plot of the frequency likelihood')
    pylab.gca().set_ylim([-100, 40])
    pylab.gca().set_autoscale_on(False)
    
for frame in stream:
    samples = frame
    fft = ffter.process( windower.process( frame ) )
    spec =  loudia.magToDb( abs( fft ) )
    
    wspec = whitening.process( spec )
    #wspec = spec
    pitch, saliency, freqs = pitchInverseProblem.process( spec )

    frameNumber -= 1
    
    if interactivePlot:
        if frameNumber == 0:
            maxFreq = freqs[0,:].argmax()
            
            pylab.subplot(211)
            pylab.hold(True)
            pylab.plot( spec[0, :plotSize] )
            pylab.plot(loudia.magToDb(a)[:, maxFreq] + 40)
            
            pylab.subplot(212)
            pylab.hold(False)
            pylab.plot(freqs[0,:plotSize], label = 'Noise Suppressed Spectrum')
            #pylab.hold(True)
            #pylab.stem( pitch/sampleRate*fftSize, saliency )
            #pylab.show()
        
    specs.append( spec[0,:plotSize] )
    wspecs.append( wspec[0,:plotSize] )
    pitches.append( pitch )
    saliencies.append( saliency )
    freqss.append( freqs )

if interactivePlot:
    pylab.ioff()


specs = scipy.array( specs )
wspecs = scipy.array( wspecs )
pitches = scipy.array( pitches )[:, 0, :]
saliencies = scipy.array( saliencies )[:, 0, :]
freqss = scipy.array( freqss )[:,0,:]
frameCount = specs.shape[0] - 1

#import scipy.signal
#print saliencies.shape
#saliencies = scipy.signal.lfilter(scipy.array([1.0/n]*int(n)), scipy.array([1.0]), saliencies, axis=0)

if plot:
    pitches[ saliencies < saliencyThreshold ] = scipy.NaN

    # Get the onsets
    annotation = os.path.splitext(filename)[0] + '.onset_annotated'
    onsets = get_onsets(annotation, frameStep, sampleRate)
    
    pylab.figure()
    pylab.subplot(311)
    pylab.imshow( scipy.flipud(freqss.T), aspect = 'auto' )

    pylab.subplot(312)
    pylab.plot( pitches[:,:] )

    draw_onsets(onsets)

    pylab.subplot(313)
    pylab.plot( saliencies[:,:] )

    draw_onsets(onsets)
    
    pylab.show()

print pitches
