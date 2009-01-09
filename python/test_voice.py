#!/usr/bin/env python

import ricaudio
from sepel.inputs import pyricaudio
import pylab
import sys
import scipy

filename = sys.argv[1]

frameSize = 1024 / 44100.0
frameStep = 512 / 44100.0

fftSize = 1024

analysisLimit = 10.0

# Creation of the pipeline        
stream = pyricaudio.sndfilereader({'filename': filename,
                                   'windowSizeInTime': frameSize,
                                   'windowStepInTime': frameStep,
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
                                             'windowType': 'blackmanharris'})

stream = pyricaudio.fft_ricaudio(stream, {'inputKey': 'windowed',
                                          'outputKey': 'fft',
                                          'zeroPhase': True,
                                          'fftLength': fftSize})



subplots = {1 : ['mag', 'peak_mags'],
            2 : ['phase', 'peak_phases']}

all_processes = set()
for values in subplots.values():
    all_processes |= set(values)

plotSize = fftSize / 4

subplotCount = max(subplots)

print subplots
pylab.ion()

if 'peak_mags' in all_processes:
    peaker = ricaudio.PeakPick( plotSize / 3 )

for frame in stream:
    fft = frame['fft'][:plotSize]
    spec =  20.0 / scipy.log(10.0) * scipy.log(abs(fft) + 1e-7)

    if set(['phase', 'peak_phases']) | all_processes:
        phase =  scipy.angle(fft)

    if set(['peak_mags', 'peak_phases']) | all_processes:
        fft = scipy.reshape(fft, (1, plotSize))
        peakLocs, peakMags =  peaker.process( fft )
        peakLocs = scipy.array(peakLocs, dtype = 'i4')

    for subplot, processes in subplots.items():
        pylab.subplot(subplotCount, 1, subplot)
        pylab.gca().clear()
        pylab.gca().set_autoscale_on(False)
        
        if 'mag' in processes:       
            pylab.gca().set_xlim([0, plotSize])
            pylab.gca().set_ylim([-100, 40])
        
            pylab.plot(spec)
            
            
        if 'peak_mags' in processes:            
            pylab.stem(peakLocs[0,:], spec[peakLocs[0,:]])
            
            
        if 'phase' in processes:           
            pylab.gca().set_xlim([0, plotSize])
            pylab.gca().set_ylim([-scipy.pi, scipy.pi])
            
            pylab.plot(phase)

            
        if 'peak_phases' in processes:
            pylab.stem(peakLocs[0,:], phase[peakLocs[0,:]])
            
pylab.ioff()
