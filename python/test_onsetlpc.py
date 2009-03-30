#!/usr/bin/env python

import loudia
from common import *
import pylab
import os, sys, wave
import scipy

interactivePlot = False

filename = sys.argv[1]

# Accepted difference between the groundtruth
# and the estimated onsets in milliseconds (ms)
onsetError = 50.0

frameSize = 1024
frameStep = 512

fftSize = 2048 * 2
plotSize = fftSize / 4

numCoeffs = 18
preEmphasis = 0.975

stream, samplerate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)

windower = loudia.Window( frameSize, loudia.Window.HAMMING )
ffter = loudia.FFT( fftSize )
lpc = loudia.LPC(frameSize, numCoeffs, preEmphasis)
lpcr = loudia.LPCResidual(frameSize)

specs = []
lpcs = []
freqResps = []
errors = []

npoints = 1024
w = scipy.arange(0, scipy.pi, scipy.pi/(fftSize/2. + 1))
b = scipy.array([[1]])

if interactivePlot:
    pylab.ion()
    pylab.figure()
    pylab.title('Interactive plot of the FFT vs LPC frequency response')
    pylab.gca().set_ylim([-100, 40])
    pylab.gca().set_autoscale_on(False)
    
for frame in stream:
    fft = ffter.process( windower.process( frame ) )[0, :]
    spec =  loudia.magToDb( abs( fft ) )
    
    lpcCoeffs, reflection, error = lpc.process( frame )

    lpcResidual = lpcr.process( frame, lpcCoeffs )
    
    freqResp = loudia.magToDb( abs( loudia.freqz( b * scipy.sqrt( abs( error[0] ) ), lpcCoeffs.T, w ) ).T )
    
    if interactivePlot:
        pylab.subplot( 211 )
        pylab.hold( False )
        pylab.plot( frame )
        #pylab.hold( True )
        #pylab.plot( lpcResidual[0,:] )
        
        pylab.subplot( 212 )
        pylab.hold( False )
        pylab.plot( w, spec[0, :], label = 'FFT' )
        pylab.hold( True )
        pylab.plot( w, freqResp[0, :], label = 'LPC' )
        
    specs.append( spec[0, :plotSize] )
    lpcs.append( lpcCoeffs[0] )
    freqResps.append( freqResp[0, :plotSize] )
    errors.append( error[0] )

if interactivePlot:
    pylab.ioff()
    
specs = scipy.array( specs )
frameCount = specs.shape[0] - 1

lpcs = scipy.array( lpcs )
freqResps = scipy.array( freqResps )
errors = scipy.array( errors )

# Get the onsets
annotation = os.path.splitext(filename)[0] + '.onset_annotated'
onsets = get_onsets(annotation, frameStep, samplerate)
    
pylab.figure()
pylab.hold(True)
pylab.subplot(3, 1, 1)

pylab.imshow( scipy.flipud(specs.T), aspect = 'auto' )

draw_onsets(onsets)
    
pylab.title( 'Spectrogram' )
ax = pylab.gca()

ax.set_xticks( ax.get_xticks()[1:] )
ticks = ax.get_xticks()
ax.set_xticklabels(['%.2f' % (float(tick) * frameStep / samplerate) for tick in ticks])

ax.set_yticks([])
ax.set_yticklabels([])

ax.set_xlim([0, frameCount - 1])


pylab.subplot(3, 1, 2)

pylab.imshow( scipy.flipud(freqResps.T), aspect = 'auto' )

pylab.title( 'LPC Frequency Response' )
ax = pylab.gca()

ax.set_xticks( ax.get_xticks()[1:] )
ticks = ax.get_xticks()
ax.set_xticklabels(['%.2f' % (float(tick) * frameStep / samplerate) for tick in ticks])

ax.set_yticks([])
ax.set_yticklabels([])

ax.set_xlim([0, frameCount - 1])


# Create the ODF processors and process
pylab.subplot(3, 1, 3)
pylab.plot( errors )

draw_onsets( onsets )

pylab.title( 'LPC Error' )
ax = pylab.gca()

ax.set_xticks([])
ax.set_yticks([])

ax.set_xticklabels([])
ax.set_yticklabels([])

ax.set_xlim([0, frameCount - 1])       

pylab.subplots_adjust(left = 0.05, right = 0.95, bottom = 0.05, top = 0.95, hspace=0.6)

pylab.show()
