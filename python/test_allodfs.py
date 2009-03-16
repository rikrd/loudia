#!/usr/bin/env python

from sepel.inputs import pyricaudio
import ricaudio
from common import *
import pylab
import os, sys, wave
import scipy

filename = sys.argv[1]

# Accepted difference between the groundtruth
# and the estimated onsets in milliseconds (ms)
onsetError = 50.0

# Samplerate of the file
wavfile = wave.open(filename,'r')
samplerate = float(wavfile.getframerate())
wavfile.close()

frameSize = 1024 
frameStep = 512

frameSizeTime = frameSize / 44100.0
frameStepTime = frameStep / 44100.0

fftSize = 2048*2
plotSize = fftSize / 8

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


odfNames = [(getattr(ricaudio.ODF, i), i) for i in dir(ricaudio.ODF)
            if type(getattr(ricaudio.ODF, i)) == int]

odfNames.sort()

specs = []
ffts = []
odfs = {}

for frame in stream:
    fft = scipy.array(frame['fft'][:fftSize/2])
    mag =  abs(fft)
    spec =  20.0 / scipy.log( 10.0 ) * scipy.log( abs( fft ) + 1e-7)[:plotSize]

    ffts.append( fft )
    specs.append( spec )

ffts = scipy.array( ffts )
specs = scipy.array( specs )

frameCount = specs.shape[0] - 1
subplots = len(odfNames) + 1

# Get the onsets
annotation = os.path.splitext(filename)[0] + '.onset_annotated'
onsets = []

if os.path.isfile(annotation):
    onsets = get_onsets(annotation, frameStep, samplerate)
        
pylab.figure()
pylab.hold(True)
pylab.subplot(subplots, 1, 1)

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
for i, (odfType, odfName) in enumerate(odfNames):
    odf = ricaudio.ODF( fftSize, odfType )

    print 'Processing: %s...' % odfName
    odfValues = odf.process( ffts )
    print 'Finished'

    pylab.subplot(subplots, 1, i+2)
    pylab.plot(odfValues[:,0])

    draw_onsets( onsets )

    pylab.title( odfName.replace('_', ' ').capitalize() )
    ax = pylab.gca()

    ax.set_xlim([0, frameCount - 1])
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_xticklabels([])
    ax.set_yticklabels([])

pylab.subplots_adjust(left = 0.05, right = 0.95, bottom = 0.05, top = 0.95, hspace=0.6)
        
pylab.show()
