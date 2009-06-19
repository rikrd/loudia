#!/usr/bin/env python

import loudia
import scipy
import pylab
import os

def plot_freqz(b, a, w = None, npoints = None, title = '', db = False, createFigure = True, label = ''):
    # Create the omega array if necessary
    if npoints is None:
        npoints = 1000

    if w is None:
        w = scipy.arange(-scipy.pi, scipy.pi, 2*scipy.pi/(npoints))

    # Calculate the frequency response
    d = loudia.freqz(b.T, a.T, w)

    if db:
        mag = 20.0 * scipy.log10(abs(d[:,0]))
    else:
        mag = abs(d[:,0])

    import pylab
    if createFigure:
        pylab.figure()
        
    pylab.subplot(2,1,1)
    pylab.plot(w, mag, label = label)
    pylab.title('%s \n Magnitude of the Frequency Response' % title)

    pylab.subplot(2,1,2)
    pylab.plot(w, scipy.angle(d[:,0]), label = label)
    pylab.title('Angle of the Frequency Response')

def get_onsets(filename, hop, sampleRate, onsetError = 50.0):
    # Get the onsets
    annotation = os.path.splitext(filename)[0] + '.onset_annotated'
    onsets = []
    
    if os.path.isfile(annotation):
        onsetsTimes = [float(o) for o in open(annotation, 'r').readlines()]
        onsetsCenter = [int(o * sampleRate / hop) for o in onsetsTimes]
        onsetsLeft = [int((o - (onsetError / 1000.0)) * sampleRate / hop) for o in onsetsTimes]
        onsetsRight = [int((o + (onsetError / 1000.0)) * sampleRate / hop) for o in onsetsTimes]
        onsets = zip(onsetsLeft, onsetsCenter, onsetsRight)

        return onsets

    else:
        return None
    
def draw_onsets(onsets):
    if not onsets:
        return
    
    # Draw the onsets
    for onsetLeft, onsetCenter, onsetRight in onsets:
        pylab.axvspan( xmin = onsetLeft, xmax = onsetRight, facecolor = 'green', linewidth = 0, alpha = 0.25)
        pylab.axvline( x = onsetCenter, color = 'black', linewidth = 1.1)


def get_framer_audio(filename, size, hop):
    from scikits import audiolab
    
    loader = audiolab.Sndfile(filename)
    sr = loader.samplerate
    nframes = loader.nframes
    nchannels = loader.channels

    framer = framer_audio(loader, size, hop)
    
    return framer, sr, nframes, nchannels, loader

def framer_audio(loader, size, hop):
    result = []
    cursor = 0L

    nchannels = loader.channels
    nframes = loader.nframes
    samples = scipy.zeros((size, nchannels))
    
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
    samples = scipy.zeros((size, arr.shape[1]))
    
    while cursor < nframes:
        nframes_read = min(size, nframes-cursor)
        samples[:nframes_read, :] = arr[cursor:cursor+nframes_read, :]
        
        # fill in empty
        if nframes_read < size:
            samples[nframes_read:, :] = 0.0

        yield samples
        cursor += hop

def overlap_add(frames, size, hop):
    nframes = len(frames)

    arrsize = size + hop * (nframes - 1)
    
    arr = scipy.zeros((arrsize, frames[0].shape[1]))

    for cur, frame in enumerate(frames):
        arr[cur*hop:(cur*hop) + size,:] += frame

    return arr
