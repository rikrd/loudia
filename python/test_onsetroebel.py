#!/usr/bin/env python

import loudia
from common import *
import pylab
import os, sys, wave
import scipy

filename = sys.argv[1]

# Accepted difference between the groundtruth
# and the estimated onsets in milliseconds (ms)
onsetError = 50.0

# Samplerate of the file
frameSize = 1024 
frameStep = 256

fftSize = 2048 * 2
plotSize = fftSize / 4

bandwidth = 4 * fftSize/frameSize
analysisLimit = scipy.inf

stream, samplerate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)

windower = loudia.Window( frameSize, loudia.Window.HAMMING )
ffter = loudia.FFT( fftSize )
odfcog = loudia.ODFCOG(fftSize, 10, bandwidth)

specs = []
cogs = []

for frame in stream:
    fft = ffter.process( windower.process( frame ) )

    spec =  loudia.magToDb( abs( fft ) )

    cog = odfcog.process( fft )
    
    specs.append( spec[0, :plotSize] )
    cogs.append( cog[0, 0] )
    

specs = scipy.array( specs )
frameCount = specs.shape[0] - 1

cogs = scipy.array( cogs )

odfRoebel = cogs

# Get the onsets
annotation = os.path.splitext(filename)[0] + '.onset_annotated'
onsets = []
if os.path.isfile(annotation):
    onsets = get_onsets(annotation, frameStep, samplerate, onsetError = onsetError)

pylab.figure()
pylab.hold(True)
pylab.subplot(2, 1, 1)

pylab.imshow( scipy.flipud(specs.T), aspect = 'auto' )

draw_onsets( onsets )
    
pylab.title( 'Spectrogram' )
ax = pylab.gca()

ax.set_xticks( ax.get_xticks()[1:] )
ticks = ax.get_xticks()
ax.set_xticklabels(['%.2f' % (float(tick) * frameStep / samplerate) for tick in ticks])

ax.set_yticks([])
ax.set_yticklabels([])

ax.set_xlim([0, frameCount - 1])


# Create the ODF processors and process
pylab.subplot(2, 1, 2)
pylab.plot( odfRoebel )

draw_onsets( onsets )

pylab.title( 'Axel Roebel''s Onset Detector' )
ax = pylab.gca()

ax.set_xticks([])
ax.set_yticks([])

ax.set_xticklabels([])
ax.set_yticklabels([])

ax.set_xlim([0, frameCount - 1])

pylab.subplots_adjust(left = 0.05, right = 0.95, bottom = 0.05, top = 0.95, hspace=0.6)
        
pylab.show()
