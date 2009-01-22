#!/usr/bin/env python

import ricaudio
from sepel.inputs import pyricaudio
import pylab
import os, sys
import scipy

filename = sys.argv[1]

samplerate = 44100

frameSize = 1024 
frameStep = 512

frameSizeTime = frameSize / 44100.0
frameStepTime = frameStep / 44100.0

fftSize = 2048
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


reassignment = ricaudio.SpectralReassignment(frameSize, fftSize, samplerate, ricaudio.Window.HAMMING)

specs = []
times = []
timesWeighted = []

for frame in stream:
    samples = scipy.array(frame['samplesMono'], dtype = 'f4')
    fft, time, freq = reassignment.process(samples)
    spec =  20.0 / scipy.log( 10.0 ) * scipy.log( abs( fft ) + 1e-7)[0, :plotSize]
    time = time[0, :plotSize]

    specs.append( spec )
    times.append( time )
    timesWeighted.append( time * spec )
    

specs = scipy.array( specs )
times = scipy.array( times )
timesWeighted = scipy.array( timesWeighted )

odfRoebel = times.sum(axis = 1) / plotSize

# Get the onsets
annotation = os.path.splitext(filename)[0] + '.onset_annotated'
onsets = []
if os.path.isfile(annotation):
    onsets = [int(float(o) * samplerate / frameStep) for o in open(annotation, 'r').readlines()]

def drawOnsets():
    # Draw the onsets
    for onset in onsets:
        pylab.axvline( x = onset, color = 'k', linewidth = 2, alpha = 0.5)

pylab.figure()
pylab.hold(True)
pylab.subplot(2, 1, 1)

pylab.imshow( scipy.flipud(specs.T), aspect = 'auto' )

drawOnsets()
    
pylab.title( 'Spectrogram' )
ax = pylab.gca()

ax.set_xticks( ax.get_xticks()[1:] )
ticks = ax.get_xticks()
ax.set_xticklabels(['%.2f' % (float(tick) * frameStep / samplerate) for tick in ticks])

ax.set_yticks([])
ax.set_yticklabels([])
    
# Create the ODF processors and process
pylab.subplot(2, 1, 2)
pylab.plot(odfRoebel)

drawOnsets()

ax = pylab.gca()

ax.set_xticks([])
ax.set_yticks([])

ax.set_xticklabels([])
ax.set_yticklabels([])

pylab.subplots_adjust(left = 0.05, right = 0.95, bottom = 0.05, top = 0.95, hspace=0.6)
        
pylab.show()
