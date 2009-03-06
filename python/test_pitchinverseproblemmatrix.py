#!/usr/bin/env python

import ricaudio
import pylab

plotFreqs = 10

pitchInverseProblem = ricaudio.PitchInverseProblem(4096, 50, 6000, 44100)
a = pitchInverseProblem.projectionMatrix()

pylab.figure()
for i in range(min(a.shape[0], plotFreqs)):
    pylab.plot(a[i,:])

pylab.figure()
for i in range(min(a.shape[1], plotFreqs)):
    pylab.plot(a[:,i])

pylab.show()
