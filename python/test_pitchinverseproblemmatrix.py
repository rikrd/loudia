#!/usr/bin/env python

import loudia
import pylab

plotFreqs = [0, 1, 2, 3, 4]
frameSize = 8192 
frameStep = 2048

fftSize = 8192

plotSize = fftSize / 8

peakBandwidth = 4
peakCandidateCount = 4
numMaxPitches = 1
numHarmonics = 10
numCandidates = 30

sampleRate = 44100

pitchInverseProblem = loudia.PitchInverseProblem(fftSize,
                                                 200.0, 4000.0,
                                                 sampleRate,
                                                 numMaxPitches,
                                                 numHarmonics,
                                                 numCandidates,
                                                 peakBandwidth)
a = pitchInverseProblem.projectionMatrix()

print "Projection matrix: ", a.shape

pylab.figure()
for i in plotFreqs:
    if i >= a.shape[0]:
        continue
    
    pylab.plot(a[i,:])
    
pylab.title("Rows: (Harmonic trains)")

pylab.figure()
for i in range(a.shape[1]):
    if i >= a.shape[1]:
        continue
    pylab.plot(a[:,i])

pylab.title("Cols: (Pitch likelihood given a bin)")

pylab.figure()
npoints = 1000
w = pylab.linspace(-100., 100., npoints)

y = pylab.zeros((1, npoints))
for i, x in enumerate(w):
    y[0, i] = loudia.gaussian(x, 0.0, 8.0)


w = pylab.array(w, dtype = 'f4')
w.resize((1, npoints))

pylab.plot(w[0,:], y[0,:])

pylab.show()
