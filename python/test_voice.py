#!/usr/bin/env python

import ricaudio
from sepel.inputs import pyricaudio
import pylab
import sys
import scipy

filename = sys.argv[1]

frameSize = 2048 / 44100.0
frameStep = 512 / 44100.0

fftSize = 2048

analysisLimit = 1000.0

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


interactivePlotting = False

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
    maxPeakCount = fftSize / 3
    maxTrajCount = fftSize / 3
    minPeakWidth = 12 # bins for Hamming
    maxFreqBinChange = 4 * fftSize / (frameSize * 44100)
    
    peaker = ricaudio.PeakDetect( maxPeakCount, minPeakWidth )
    peakInterp = ricaudio.PeakInterpolate( )
    tracker = ricaudio.PeakContinue( maxTrajCount, maxFreqBinChange )

trajsLocs = []
trajsMags = []
specs = []

for frame in stream:
    fft = scipy.array(frame['fft'][:plotSize], dtype = scipy.complex64)
    mag =  scipy.array(abs(fft), dtype = 'f4')
    spec =  20.0 / scipy.log( 10.0 ) * scipy.log( abs( fft ) + 1e-7)

    if set(['phase', 'peak_phases']) | all_processes:
        phase =  scipy.angle( fft )

    if set(['peak_mags', 'peak_phases']) | all_processes:
        fft = scipy.reshape(fft, (1, plotSize))
        peakLocs, peakMags =  peaker.process( fft )
        peakiLocs, peakiMags = peakInterp.process( fft,
                                                   scipy.array(peakLocs, dtype='f4'),
                                                   scipy.array(peakMags, dtype='f4'))
        
        trajLocs, trajMags = tracker.process( fft,
                                              scipy.array(peakiLocs, dtype='f4'),
                                              scipy.array(peakiMags, dtype='f4') )

        
        trajsLocs.append( trajLocs[0,:] )
        trajsMags.append( trajMags[0,:] )

        specs.append( spec )

        peakPos = peakLocs[peakLocs > 0]
        peakMags = peakMags[peakLocs > 0]

    if interactivePlotting:
        for subplot, processes in subplots.items():
            pylab.subplot(subplotCount, 1, subplot)
            pylab.gca().clear()
            pylab.gca().set_autoscale_on(False)

            if 'mag' in processes:       
                pylab.gca().set_xlim([0, plotSize])
                pylab.gca().set_ylim([-100, 40])

                pylab.plot(spec)


            if 'peak_mags' in processes:
                pylab.scatter(peakPos, spec[scipy.array(peakPos, dtype='i4')], c='r')


            if 'phase' in processes:           
                pylab.gca().set_xlim([0, plotSize])
                pylab.gca().set_ylim([-scipy.pi, scipy.pi])

                pylab.plot(phase)


            if 'peak_phases' in processes:
                pylab.scatter(peakPos, phase[scipy.array(peakPos, dtype='i4')], c='r')
    
            
            
pylab.ioff()

"""
trajsLocs = scipy.array( trajsLocs )
trajsMags = scipy.array( trajsMags )


def extractTrajs(trajsLocs, trajsMags):
    trajs = []
    
    for col in range(trajsLocs.shape[1]):
        inTrack = False
        trackInds = []
        trackMags = []
        trackPos = []
        
        for row in range(trajsLocs.shape[0]):
            if not trajsMags[row, col] == 0:
                inTrack = True
                
                trackInds.append( row )
                trackMags.append( trajsMags[row, col] )
                trackPos.append( trajsLocs[row, col] )
            else:
                if inTrack:
                    
                    trackInds = scipy.array(trackMags)
                    trackMags = scipy.array(trackMags)
                    trackPos = scipy.array(trackPos)
                    trajs.append((trackInds, trackPos, trackMags))

                inTrack = False
                    
                trackInds = []
                trackMags = []
                trackPos = []


        if inTrack:
            trackInds = scipy.array(trackMags)
            trackMags = scipy.array(trackMags)
            trackPos = scipy.array(trackPos)
            trajs.append((trackInds, trackPos, trackMags))

        

    print len(trajs)
    return trajs


trajs = extractTrajs(trajsLocs, trajsMags)

pylab.figure(2)
for trajInds, trajPos, trajMags in trajs:
    pylab.hold( True )
    pylab.scatter( trajInds, trajPos )
"""
trajsLocs = scipy.array( trajsLocs )
trajsMags = scipy.array( trajsMags )

print trajsLocs.shape

pylab.figure()
pylab.plot( trajsLocs )

pylab.figure()
pylab.plot( trajsMags )

pylab.figure()
pylab.hold(True)
pylab.imshow( scipy.array( specs ).T, aspect = 'equal' )
pylab.plot( trajsLocs )

pylab.show()
