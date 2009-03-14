#!/usr/bin/env python

import ricaudio
from common import *
import pylab
import scipy.signal

samplerate = 22050.0

order = 15
center = 1000.0
bandwidth = 200.0
ripplePass = 0.05
rippleStop = 40
npoints = 200

btype = 'bandstop'
ftype = 'butter'

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
for center in range(2000.0, 5000.0, 1000.0):
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
    
    rb = f.b().T#[:,:1]
    ra = f.a().T

    plotFreqz(rb, ra,
              npoints = npoints,
              createFigure = False,
              label = 'ricaudio / %.1f Hz' % bandwidth,
              db = True)

    # For Scipy
    sb, sa = scipy.signal.iirfilter(order, [freq, freqStop],
                                    btype = btype, ftype = ftype,
                                    rp = rippleStop, rs = ripplePass)
    
    sa = scipy.array(scipy.reshape(sa, (1, sa.shape[0])), dtype = 'f4')
    sb = scipy.array(scipy.reshape(sb, (1, sb.shape[0])), dtype = 'f4')#[:,:1]

    results[bandwidth] = scipy.allclose(sa, ra) and scipy.allclose(sb, rb)
    print ra
    print sa
    print rb
    print sb
    #diffs[bandwidth] = (ra - sa, rb - sb)
    
    plotFreqz(sb, sa,
              npoints = npoints,
              createFigure = False,
              label = 'scipy / %.1f Hz' % bandwidth,
              db = True)

print str(results)
print str(diffs)

pylab.legend()
pylab.show()
