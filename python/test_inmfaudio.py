#!/usr/bin/env python

from common import *
import sys

plot = True
plotInteractive = False

estimatePitch = True

filename = '/home/rmarxer/dev/data/onsets/pitchedphrases/Strings/Picked-Plucked-Ham/Piano/piano1.wav'

if len(sys.argv) >= 2:
    filename = sys.argv[1]

windowSize = 1024
windowHop = 256
windowType = ricaudio.Window.BLACKMANHARRIS
fftSize = 1024
halfSize = (fftSize / 2) + 1
plotSize = halfSize/3
zeroPhase = True

components = 5
pastFrames = 20
pastCoeff = 0.2

w = ricaudio.Window(windowSize, windowType)
f = ricaudio.FFT(windowSize, fftSize, zeroPhase)
d = ricaudio.INMF(halfSize, components, pastFrames, pastCoeff, 30, 1, 1e-17)

framer, sr, nframes, nchannels, loader = get_framer_audio(filename, windowSize, windowHop)

onsets = get_onsets(filename, windowHop, sr)

if estimatePitch:
   acorrw = ricaudio.Window(halfSize, ricaudio.Window.HAMMING)
   peaker = ricaudio.PeakDetect(1, 3, 0, False)
   peakeri = ricaudio.PeakInterpolate()
   freqs = []

components = []
gains = []
a = []
for frame in framer:
    windowed = w.process(frame)
    fft = f.process(windowed)
    fft = abs(fft)
    c, g = d.process(fft)
    
    fft = 20.0 * scipy.log10(fft[0,:plotSize] + 0.1)
    dBc = 20.0 * scipy.log10(c[:,:plotSize].T + 0.1)

    if estimatePitch:
        acorrc = acorrw.process( c )
        acorr = ricaudio.autocorrelate( acorrc )
        peakPos, peakMag, peakPhase = peaker.process( acorr )
        peakPosi, peakMagi, peakPhasei = peakeri.process( acorr, peakPos, peakMag, peakPhase )
        
        fs = peakPosi / fftSize * sr

        if plotInteractive:
            pylab.ion()
            pylab.figure( 1 )
            pylab.subplot( 211 )
            pylab.hold( False )
            pylab.plot( acorrc.T )
            
            pylab.subplot( 212 )
            pylab.plot( acorr.T )

        freqs.append( fs.T )
            
    if plotInteractive:
        pylab.ion()
        pylab.figure( 2 )
        pylab.hold( False )
        pylab.plot( dBc )
        
    a.append( fft )
    gains.append( g )
    components.append( dBc )

pylab.ioff()

gains = overlapadder(gains, 1, 1)
lastComponents = components[-1]
a = scipy.array(a)
nwindows = a.shape[0]

if estimatePitch:
    freqs = overlapadder(freqs, 1, 1)

    dBc = scipy.array(components)
    dBcchange = abs(dBc[2:,:,:] - dBc[1:-1,:,:])
    odfComponents = (gains[2:,:] * dBcchange.sum(axis = 1)).sum(axis = 1)
    
    freqchange = abs(freqs[2:,:] - freqs[1:-1,:])
    wfreqchange = gains[2:,:] * freqchange
    odfPitch = wfreqchange.sum(axis = 1) #/ gains[1:, :].sum(axis = 1)

    if plot:
        pylab.figure()
        pylab.subplot( 311 )
        pylab.plot( freqs )
        draw_onsets( onsets )        
        pylab.title("Pitch Contours")
        pylab.gca().set_xlim([0, nwindows - 1])

        pylab.subplot( 312 )
        pylab.plot( odfPitch )
        draw_onsets( onsets )        
        pylab.title("INMF Pitch derived ODF")
        pylab.gca().set_xlim([0, nwindows - 1])

        pylab.subplot( 313 )
        pylab.plot( odfComponents )
        draw_onsets( onsets )        
        pylab.title("INMF Components derived ODF")
        pylab.gca().set_xlim([0, nwindows - 1])


if plot:
    pylab.figure()
    pylab.subplot(211)
    pylab.imshow(scipy.flipud(a.T), aspect = 'auto')
    draw_onsets(onsets)
    pylab.title("Input")
    pylab.gca().set_xlim([0, nwindows - 1])


    pylab.subplot(212)
    pylab.plot(gains)
    draw_onsets(onsets)
    pylab.title("Gains")
    pylab.gca().set_xlim([0, nwindows - 1])


    pylab.figure()
    pylab.plot(lastComponents)
    pylab.title("Last Components")

    pylab.show()
