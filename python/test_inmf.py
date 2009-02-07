#!/usr/bin/env python

import ricaudio
import scipy
import pylab

plot = True

windowSize = 14
windowHop = 14

d = ricaudio.INMF(8, 2, 10, 0.5, 15, 1, 1e-17)


def framer(arr, size, hop):
    result = []
    cursor = 0L

    nframes = arr.shape[0]
    samples = scipy.zeros((size, arr.shape[1]), dtype = a.dtype)
    
    while cursor < nframes:
        nframes_read = min(size, nframes-cursor)
        samples[:nframes_read, :] = arr[:nframes_read, :]
        
        # fill in empty
        if nframes_read < size:
            samples[nframes_read:, :] = 0.0

        yield samples
        cursor += hop

def overlapadder(frames, size, hop):
    nframes = len(frames)

    arrsize = size + hop * (nframes - 1)

    arr = scipy.zeros((arrsize, frames[0].shape[1]))

    for cur, frame in enumerate(frames):
        print frame.shape
        print cur*hop + size
        arr[cur*hop:(cur*hop) + size,:] += frame

    return arr

a = scipy.zeros((14, 8), dtype = 'f4')
a[:4, 2] = 1
a[5:9, 6] = 1
a[10:, 2] = 1

components = []
gains = []
for frame in framer(a, windowSize, windowHop):
    c, g = d.process(frame)

    gains.append(g)
    components.append(c)

gains = overlapadder(gains, windowSize, windowHop)
components = components[-1]

if plot:
    pylab.figure()
    pylab.imshow(scipy.flipud(a.T))
    pylab.title("Input")
    
    pylab.figure()
    pylab.plot(components.T)
    pylab.title("Components")
    
    pylab.figure()
    pylab.plot(gains)
    pylab.title("Gains")
    
    pylab.show()
