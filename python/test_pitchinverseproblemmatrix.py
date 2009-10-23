#!/usr/bin/env python

import loudia
import pylab

frameSize = 8192 
frameStep = 2048

fftSize = 8192

plotSize = fftSize / 8

peakBandwidth = 3
peakCandidateCount = 4
numMaxPitches = 1
numHarmonics = 10
numCandidates = 300

sampleRate = 44100

pitchInverseProblem = loudia.PitchInverseProblem(fftSize,
                                                 200.0, 4000.0,
                                                 sampleRate,
                                                 numMaxPitches,
                                                 numHarmonics,
                                                 numCandidates,
                                                 peakBandwidth)
a = pitchInverseProblem.projectionMatrix()

nBins = 10
plotBins = map(lambda x: x/float(nBins)*a.shape[0], range(nBins))

nFreqs = 4
plotFreqs = map(lambda x: x/float(nFreqs)*a.shape[1], range(nFreqs))

print "Projection matrix: ", a.shape

pylab.figure()
for i in plotBins:
    if i >= a.shape[0]:
        continue
    
    pylab.plot(a[i,:], label="bin: %d" % i)
    
pylab.title("Rows: (Pitch likelihood given a bin)")
pylab.legend()

pylab.figure()
for i in plotFreqs:
    if i >= a.shape[1]:
        continue
    pylab.plot(a[:,i], label="freq: %d" % i)

pylab.title("Cols: (Harmonic trains)")
pylab.legend()

## pylab.figure()
## npoints = 1000
## w = pylab.linspace(-100., 100., npoints)

## y = pylab.zeros((1, npoints))
## for i, x in enumerate(w):
##     y[0, i] = loudia.gaussian(x, 0.0, 8.0)


## w = pylab.array(w, dtype = 'f4')
## w.resize((1, npoints))

## pylab.plot(w[0,:], y[0,:])

pylab.show()
