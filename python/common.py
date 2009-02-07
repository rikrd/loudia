#!/usr/bin/env python

import ricaudio
import scipy
import pylab


def framer(arr, size, hop):
    result = []
    cursor = 0L

    nframes = arr.shape[0]
    samples = scipy.zeros((size, arr.shape[1]), dtype = arr.dtype)
    
    while cursor < nframes:
        nframes_read = min(size, nframes-cursor)
        samples[:nframes_read, :] = arr[cursor:cursor+nframes_read, :]
        
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
        arr[cur*hop:(cur*hop) + size,:] += frame

    return arr
