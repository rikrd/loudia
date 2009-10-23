#!/usr/bin/env python

import loudia
from common import *
import pylab
import os, sys, wave
import scipy


interactivePlot = False
plot = True

filename = sys.argv[1]

frameSize = 8192 
frameStep = 2048

fftSize = 8192

plotSize = fftSize / 8

stream, sampleRate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)


peakBandwidth = 3
peakCandidateCount = 4
numMaxPitches = 1
numHarmonics = 80
numCandidates = 300

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

specs = []
wspecs = []
pitches = []
saliencies = []
freqss = []

if interactivePlot:
    pylab.ion()
    pylab.figure()
    pylab.title('Interactive plot of the FFT vs LPC frequency response')
    pylab.gca().set_ylim([-100, 40])
    pylab.gca().set_autoscale_on(False)

for frame in stream:
    samples = frame
    fft = ffter.process( windower.process( frame ) )
    spec =  loudia.magToDb( abs( fft ) )
    
    wspec = whitening.process( spec )
    #wspec = spec
    pitch, saliency, freqs = pitchInverseProblem.process( wspec )
    
    if interactivePlot:
        pylab.subplot(211)
        pylab.hold(False)
        pylab.plot( wspec[0, :plotSize] )

        pylab.subplot(212)
        pylab.hold(False)
        pylab.plot(freqs[0,:plotSize], label = 'Noise Suppressed Spectrum')
        #pylab.hold(True)
        #pylab.stem( pitch/sampleRate*fftSize, saliency )
        
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
    pitches[ saliencies < 35] = scipy.NaN

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
