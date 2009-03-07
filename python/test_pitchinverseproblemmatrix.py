#!/usr/bin/env python

import ricaudio
import pylab

plotFreqs = 100

pitchInverseProblem = ricaudio.PitchInverseProblem(4096, 100, 5100, 44100, 5, 10, 512, 4)
a = pitchInverseProblem.projectionMatrix()

pylab.figure()
for i in range(min(a.shape[0], plotFreqs)):
    pylab.plot(a[i,:])

pylab.figure()
for i in range(min(a.shape[1], plotFreqs)):
    pylab.plot(a[:,i])


pylab.figure()
npoints = 1000
w = pylab.linspace(-100., 100., npoints)

y = pylab.zeros((1, npoints))
for i, x in enumerate(w):
    y[0, i] = ricaudio.gaussian(x, 0.0, 8.0)


w = pylab.array(w, dtype = 'f4')
w.resize((1, npoints))

pylab.plot(w[0,:], y[0,:])

pylab.show()
