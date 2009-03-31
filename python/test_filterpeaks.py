#!/usr/bin/env python

import loudia
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

order = 10

filename = sys.argv[1]

frameSize = 1024 
frameStep = 256

fftSize = 2048

stream, samplerate, nframes, nchannels, loader = get_framer_audio(filename, frameSize, frameStep)

ffter = loudia.FFT( fftSize )
windower = loudia.Window( frameSize, loudia.Window.HAMMING )

subplots = {1 : ['mag', 'peaki_mags', 'resid_mag', 'synth_mag', 'traj_mags'],
            2 : ['phase', 'peak_phases']}

all_processes = set()
for values in subplots.values():
    all_processes |= set(values)

plotSize = fftSize / 4

subplotCount = max(subplots)

pylab.ion()

if 'peaki_mags' in all_processes:
    maxPeakCount = 30
    maxTrajCount = 30
    silentFrames = 3
    minPeakWidth = 4 * int(fftSize / frameSize) # bins for Hamming
    minPeakContrast = 0.0
    maxFreqBinChange = 1 * fftSize / frameSize
    windowType = loudia.Window.HAMMING
    
    peaker = loudia.PeakDetectionComplex( maxPeakCount, loudia.PeakDetectionComplex.BYMAGNITUDE, minPeakWidth )
    peakInterp = loudia.PeakInterpolationComplex( )
    tracker = loudia.PeakTracking( maxTrajCount, maxFreqBinChange, silentFrames )
    filt = loudia.BandFilter( )

trajsLocs = []
trajsMags = []
specs = []
filtereds = []

for frame in stream:
    samples = frame
    fft = ffter.process( windower.process( frame ) )[0, :]
    spec =  loudia.magToDb( abs( fft ) )[0, :plotSize]

    if set(['phase', 'peak_phases']) | all_processes:
        phase =  scipy.angle( fft )

    if set(['peak_mags', 'peak_phases']) | all_processes:
        
        peakLocs, peakMags, peakPhases =  peaker.process( fft )

        peakiLocs, peakiMags, peakiPhases = peakInterp.process( fft,
                                                                peakLocs,
                                                                peakMags,
                                                                peakPhases )
        
        trajLocs, trajMags = tracker.process( fft,
                                              peakiLocs,
                                              peakiMags )
        

        filtered = scipy.reshape(samples, (samples.shape[0], 1))

        for trajLoc in trajLocs[0, :]:
            if trajLoc == -1:
                continue
            
            # Select the critical frequencies for the given trajectory
            freq = trajLoc / fftSize - 0.01
            freqStop = trajLoc / fftSize + 0.01

            # Create the filter for the given trajectory
            filt.setOrder( order, False )
            filt.setLowFrequency( freq, False )
            filt.setHighFrequency( freqStop, False )
            filt.setFilterType( loudia.BandFilter.BESSEL, False )
            filt.setBandType( loudia.BandFilter.BANDSTOP, False )
            filt.setup()
            
            # Filter the samples of that trajectory
            filtered = filt.process( filtered )
            
        filtereds.append( filtered.T )
        
        
        specMag = scipy.resize(spec, (1, spec.shape[0]))
        
        trajsLocs.append( trajLocs[0,:] )
        trajsMags.append( trajMags[0,:] )

        specs.append( spec )

        peakPos = peakLocs[peakLocs > 0]
        peakMags = peakMags[peakLocs > 0]

        peakiPos = peakiLocs[peakiLocs > 0]
        peakiMags = peakiMags[peakiLocs > 0]

        trajPos = trajLocs[trajLocs > 0]
        trajMags = trajMags[trajLocs > 0]

    if interactivePlotting:
        for subplot, processes in subplots.items():
            pylab.subplot(subplotCount, 1, subplot)
            pylab.gca().clear()
            pylab.gca().set_autoscale_on(False)

            if 'mag' in processes:       
                pylab.gca().set_xlim([0, plotSize])
                pylab.gca().set_ylim([-100, 40])

                pylab.plot(spec)

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

filtereds = scipy.array( filtereds )

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

specs = scipy.array( specs )

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
        
pylab.show()
