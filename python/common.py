#!/usr/bin/env python

import ricaudio
import scipy
import pylab

def get_framer_audio(filename, size, hop):
    from scikits import audiolab
    
    loader = audiolab.sndfile(filename)
    sr = loader.get_samplerate()
    nframes = loader.get_nframes()
    nchannels = loader.get_channels()

    framer = framer_audio(loader, size, hop)
    
    return framer, sr, nframes, nchannels, loader

def framer_audio(loader, size, hop):
    result = []
    cursor = 0L

    nchannels = loader.get_channels()
    nframes = loader.get_nframes()
    samples = scipy.zeros((size, nchannels), dtype = 'f4')
    
    while cursor < nframes:
        nframes_read = min(size, nframes-cursor)

        loader.seek(cursor)

        if nchannels == 1:
            samples[:nframes_read, 0] = loader.read_frames(nframes_read)
        else:
            samples[:nframes_read, :] = loader.read_frames(nframes_read)[::]
        
        # fill in empty
        if nframes_read < size:
            samples[nframes_read:, :] = 0.0
        
        yield samples.mean(axis = 1).T
        cursor += hop

def framer_array(arr, size, hop):
    result = []
    cursor = 0L

    nframes = arr.shape[0]
    samples = scipy.zeros((size, arr.shape[1]), dtype = 'f4')
    
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
