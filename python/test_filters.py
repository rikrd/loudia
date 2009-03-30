#!/usr/bin/env python

import loudia
from common import *
import pylab
import scipy.signal

samplerate = 44100.0

order = 15
center = 1000.0
bandwidth = 200.0
ripplePass = 0.05
rippleStop = 40
npoints = 1000

btype = 'bandstop'
ftype = 'butter'

btypes = {'lowpass': loudia.LowPass,
          'highpass': loudia.HighPass,
          'bandpass': loudia.BandPass,
          'bandstop': loudia.BandStop}

ftypes = {'cheby1': loudia.CHEBYSHEVI,
          'cheby2': loudia.CHEBYSHEVII,
          'butter': loudia.BUTTERWORTH,
          'bessel': loudia.BESSEL}

results = {}
diffs = {}

pylab.figure()
pylab.hold( True )
for center in range(2000.0, 5000.0, 1000.0):
    freq = (center - bandwidth / 2.0) / samplerate
    freqStop = (center + bandwidth / 2.0) / samplerate

    # For Ricaudio
    loudiaFilter = btypes[btype]

    if btype.startswith('band'):
        f = loudiaFilter(order, freq, freqStop, ftypes[ftype],
                           ripplePass, rippleStop)
    else:
        f = loudiaFilter(order, freq, ftypes[ftype],
                           ripplePass, rippleStop)
    
    rb = f.b().T#[:,:1]
    ra = f.a().T

    plot_freqz(rb, ra,
               npoints = npoints,
               createFigure = False,
               label = 'loudia / %.1f Hz' % bandwidth,
               db = True)

    # For Scipy
    sb, sa = scipy.signal.iirfilter(order, [freq, freqStop],
                                    btype = btype, ftype = ftype,
                                    rp = ripplePass, rs = rippleStop)

    
    results[bandwidth] = scipy.allclose(sa, ra) and scipy.allclose(sb, rb)
    print ra
    print sa
    print rb
    print sb
    #diffs[bandwidth] = (ra - sa, rb - sb)
    
    plot_freqz(sb, sa,
               npoints = npoints,
               createFigure = False,
               label = 'scipy / %.1f Hz' % bandwidth,
               db = True)

print str(results)
print str(diffs)

pylab.legend()
pylab.show()
