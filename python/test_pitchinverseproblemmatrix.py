#!/usr/bin/env python

import ricaudio
import pylab

pitchInverseProblem = ricaudio.PitchInverseProblem(4096, 50, 6000, 44100)
a = pitchInverseProblem.projectionMatrix()

pylab.figure()
for i in range(min(a.shape[0], 20)):
    pylab.plot(a[i,:])

pylab.show()
