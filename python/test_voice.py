#!/usr/bin/env python

import ricaudio
from common import *
import pylab
import sys
import scipy


interactivePlotting = False
plotSpectrumTrajs = True
plotDetSpecSynth = True
plotDetSpecDiff = True
plotTrajs = True
plotMags = False
plotLocs = False
plotOdf = True

filename = sys.argv[1]

frameSize = 1024 
frameStep = 256

fftSize = 2048

stream, samplerate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)

windower = ricaudio.Window( frameSize, ricaudio.Window.HAMMING )
ffter = ricaudio.FFT( fftSize )

subplots = {1 : ['mag', 'peaki_mags', 'resid_mag', 'synth_mag', 'traj_mags'],
            2 : ['phase', 'peak_phases']}

all_processes = set()
for values in subplots.values():
    all_processes |= set(values)

plotSize = fftSize / 4

subplotCount = max(subplots)

print subplots
pylab.ion()

if 'peaki_mags' in all_processes:
    maxPeakCount = 10
    maxTrajCount = 20
    silentFrames = 3
    minPeakWidth = 4 * int(fftSize / frameSize) # bins for Hamming
    minPeakContrast = 0.0
    maxFreqBinChange = 1 * fftSize / frameSize
    windowType = ricaudio.Window.HAMMING
    
    peaker = ricaudio.PeakDetectionComplex( maxPeakCount, ricaudio.PeakDetectionComplex.BYMAGNITUDE, minPeakWidth )
    peakInterp = ricaudio.PeakInterpolationComplex( )
    tracker = ricaudio.PeakTracking( maxTrajCount, maxFreqBinChange, silentFrames )
    peakSynth = ricaudio.PeakSynthesize( frameSize/6, fftSize, windowType )

trajsLocs = []
trajsMags = []
specs = []
specsSynth = []
specsResid = []
specsMagsResid = []

for frame in stream:
    samples = frame
    fft = ffter.process( windower.process( frame ) )[0, :plotSize]
    spec =  ricaudio.magToDb( abs( fft ) )[0, :plotSize]

    if set(['phase', 'peak_phases']) | all_processes:
        phase =  scipy.angle( fft )

    if set(['peak_mags', 'peak_phases']) | all_processes:
        fft = scipy.reshape(fft, (1, plotSize))

        peakLocs, peakMags, peakPhases =  peaker.process( fft )

        peakiLocs, peakiMags, peakiPhases = peakInterp.process( fft,
                                                               peakLocs,
                                                               peakMags,
                                                               peakPhases )
        
        trajLocs, trajMags = tracker.process( fft,
                                              peakiLocs,
                                              peakiMags )

        specSynth = peakSynth.process( trajLocs,
                                       trajMags )

        specSynth = specSynth[:,:plotSize]

        specMag = scipy.resize(spec, (1, spec.shape[0]))

        specMagResid = ricaudio.dbToMag( specMag ) - specSynth
        
        specResid = ricaudio.magToDb(ricaudio.dbToMag( specMag ) - specSynth)[0,:]
        
        specSynth = ricaudio.magToDb( specSynth )[0,:]
        
        trajsLocs.append( trajLocs[0,:] )
        trajsMags.append( trajMags[0,:] )

        specs.append( spec )

        specsSynth.append( specSynth )
        specsResid.append( specResid )

        peakPos = peakLocs[peakLocs > 0]
        peakMags = peakMags[peakLocs > 0]

        peakiPos = peakiLocs[peakiLocs > 0]
        peakiMags = peakiMags[peakiLocs > 0]

        trajPos = trajLocs[trajLocs > 0]
        trajMags = trajMags[trajLocs > 0]

        specsMagsResid.append( specMagResid[0,:] )

    if interactivePlotting:
        for subplot, processes in subplots.items():
            pylab.subplot(subplotCount, 1, subplot)
            pylab.gca().clear()
            pylab.gca().set_autoscale_on(False)

            if 'mag' in processes:       
                pylab.gca().set_xlim([0, plotSize])
                pylab.gca().set_ylim([-100, 40])

                pylab.plot(spec)

            if 'synth_mag' in processes:       
                pylab.gca().set_xlim([0, plotSize])
                pylab.gca().set_ylim([-100, 40])

                pylab.plot(specSynth)

            if 'resid_mag' in processes:       
                pylab.gca().set_xlim([0, plotSize])
                pylab.gca().set_ylim([-100, 40])

                pylab.plot(specResid)


            if 'peaki_mags' in processes:
                if not (peakiPos == -1).all():
                    pylab.scatter(peakiPos, peakiMags, c='r')

            if 'traj_mags' in processes:
                if not (trajPos == -1).all():
                    pylab.scatter(trajPos, trajMags, c='g')


            if 'phase' in processes:           
                pylab.gca().set_xlim([0, plotSize])
                pylab.gca().set_ylim([-scipy.pi, scipy.pi])

                pylab.plot(phase)


            if 'peak_phases' in processes:
                if not (peakPos == -1).all():
                    pylab.scatter(peakPos, phase[peakPos], c='r')
    
            
            
pylab.ioff()


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
            if not trajsMags[row, col] == -120:
                inTrack = True
                
                trackInds.append( row )
                trackMags.append( trajsMags[row, col] )
                trackPos.append( trajsLocs[row, col] )
            else:
                if inTrack:
                    
                    trackInds = scipy.array(trackInds)
                    trackMags = scipy.array(trackMags)
                    trackPos = scipy.array(trackPos)
                    trajs.append((trackInds, trackPos, trackMags))

                inTrack = False
                    
                trackInds = []
                trackMags = []
                trackPos = []


        if inTrack:
            trackInds = scipy.array(trackInds)
            trackMags = scipy.array(trackMags)
            trackPos = scipy.array(trackPos)
            trajs.append((trackInds, trackPos, trackMags))

        

    print 'Number of trajectories:', len(trajs)
    return trajs


trajs = extractTrajs(trajsLocs, trajsMags)


specsSynth = scipy.array( specsSynth )
specs = scipy.array( specs )
specsDiff = scipy.array( specsResid )
specsMagsResid = scipy.array( specsMagsResid )

if plotDetSpecSynth:
    pylab.figure()
    pylab.imshow( specsSynth.T )

if plotDetSpecDiff:
    pylab.figure()
    pylab.imshow( specsDiff.T )

if plotTrajs:
    pylab.figure()
    for trajInds, trajPos, trajMags in trajs:
        pylab.hold( True )
        pylab.plot( trajInds, trajPos )

if plotLocs:
    pylab.figure()
    trajsLocs = scipy.array( trajsLocs )
    pylab.plot( trajsLocs )

if plotMags:
    pylab.figure()
    trajsMags = scipy.array( trajsMags )
    trajsMags = trajsMags.sum(axis = 1)
    pylab.plot( trajsMags )

if plotSpectrumTrajs:
    pylab.figure()
    pylab.hold(True)
    pylab.imshow( specs.T, aspect = 'auto' )
    for trajInds, trajPos, trajMags in trajs:
        pylab.hold( True )
        pylab.plot( trajInds, trajPos, c='black' )

if plotOdf:
    odf = ricaudio.ODFComplex( fftSize )
    odfValue = []
    for i in range(specsMagsResid.shape[0] - 10):
        val = odf.process(specsMagsResid[i:i+3,:])[0,0]
        odfValue.append(val)

    pylab.figure()
    pylab.plot(odfValue)
        
pylab.show()
