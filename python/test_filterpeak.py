#!/usr/bin/env python

import ricaudio
from common import *
import pylab
import scipy.signal

samplerate = 44100.0

order = 15
center = 400.0
bandwidth = 20.0
ripplePass = 0.05
rippleStop = 40

btype = 'bandstop'
ftype = 'cheby2'

btypes = {'lowpass': ricaudio.LowPass,
          'highpass': ricaudio.HighPass,
          'bandpass': ricaudio.BandPass,
          'bandstop': ricaudio.BandStop}

ftypes = {'cheby1': ricaudio.CHEBYSHEVI,
          'cheby2': ricaudio.CHEBYSHEVII,
          'butter': ricaudio.BUTTERWORTH,
          'bessel': ricaudio.BESSEL}

results = {}
diffs = {}

pylab.figure()
pylab.hold( True )
for bandwidth in range(20.0, 50.0, 10.0):
    freq = (center - bandwidth / 2.0) / samplerate
    freqStop = (center + bandwidth / 2.0) / samplerate

    # For Ricaudio
    ricaudioFilter = btypes[btype]

    if btype.startswith('band'):
        f = ricaudioFilter(order, freq, freqStop, ftypes[ftype],
                           ripplePass, rippleStop)
    else:
        f = ricaudioFilter(order, freq, ftypes[ftype],
                           ripplePass, rippleStop)
    
    rb = f.b().T[:,:1]
    ra = f.a().T

    plotFreqz(rb, ra,
              createFigure = False,
              label = 'ricaudio / %.1f Hz' % bandwidth,
              db = True)

    # For Scipy
    sb, sa = scipy.signal.iirfilter(order, [freq, freqStop],
                                    btype = btype, ftype = ftype,
                                    rp = rippleStop, rs = ripplePass)
    
    sa = scipy.array(scipy.reshape(sa, (1, sa.shape[0])), dtype = 'f4')
    sb = scipy.array(scipy.reshape(sb, (1, sb.shape[0])), dtype = 'f4')[:,:1]

    results[bandwidth] = scipy.allclose(sa, ra) and scipy.allclose(sb, rb)
    print ra
    print sa
    print rb
    print sb
    #diffs[bandwidth] = (ra - sa, rb - sb)
    
    plotFreqz(sb, sa,
              createFigure = False,
              label = 'scipy / %.1f Hz' % bandwidth,
              db = True)

print str(results)
print str(diffs)

pylab.legend()
pylab.show()
