#!/usr/bin/env python

import ricaudio
import scipy
import pylab
import os

def get_onsets(filename, hop, samplerate, onsetError = 50.0):
    # Get the onsets
    annotation = os.path.splitext(filename)[0] + '.onset_annotated'
    onsets = []
    
    if os.path.isfile(annotation):
        onsetsTimes = [float(o) for o in open(annotation, 'r').readlines()]
        onsetsCenter = [int(o * samplerate / hop) for o in onsetsTimes]
        onsetsLeft = [int((o - (onsetError / 1000.0)) * samplerate / hop) for o in onsetsTimes]
        onsetsRight = [int((o + (onsetError / 1000.0)) * samplerate / hop) for o in onsetsTimes]
        onsets = zip(onsetsLeft, onsetsCenter, onsetsRight)

        return onsets

    else:
        return None
    
def draw_onsets(onsets):
    # Draw the onsets
    for onsetLeft, onsetCenter, onsetRight in onsets:
        pylab.axvspan( xmin = onsetLeft, xmax = onsetRight, facecolor = 'green', linewidth = 0, alpha = 0.25)
        pylab.axvline( x = onsetCenter, color = 'black', linewidth = 1.1)


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
