#!/usr/bin/env python

import ricaudio
import scipy
import scipy.signal
from common import *

plot = True

atol = 1e-4

freq = 0.2
freqStop = 0.4
fs = 8000
rp = 0.05
rs = 40

### Chebyshev I
# Test with even order
order = 4

rc = ricaudio.BandPass( order, freq, freqStop, ricaudio.CHEBYSHEVI, rp, rs )
sc_b, sc_a = scipy.signal.cheby1(order, rp, [freq, freqStop], btype='bandpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'BandPass Chebyshev I order %d' % order)

# Test with odd order
order = 5

rc = ricaudio.BandPass( order, freq, freqStop, ricaudio.CHEBYSHEVI, rp, rs )
sc_b, sc_a = scipy.signal.cheby1( order, rp, [freq, freqStop], btype='bandpass', analog=0, output='ba' )

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'BandPass Chebyshev I order %d' % order)


### Chebyshev II
# Test with even order
order = 4

rc = ricaudio.BandPass( order, freq, freqStop, ricaudio.CHEBYSHEVII, rp, rs )
sc_b, sc_a = scipy.signal.cheby2(order, rs, [freq, freqStop], btype='bandpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'BandPass Chebyshev II order %d' % order)


# Test with odd order
order = 5

rc = ricaudio.BandPass( order, freq, freqStop, ricaudio.CHEBYSHEVII, rp, rs )
sc_b, sc_a = scipy.signal.cheby2(order, rs, [freq, freqStop], btype='bandpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'BandPass Chebyshev II order %d' % order)


### Butterworth
# Test with even order
order = 4

rc = ricaudio.BandPass( order, freq, freqStop, ricaudio.BUTTERWORTH, rp, rs )
sc_b, sc_a = scipy.signal.butter(order, [freq, freqStop], btype='bandpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'BandPass Butterworth order %d' % order)


# Test with odd order
order = 5

rc = ricaudio.BandPass( order, freq, freqStop, ricaudio.BUTTERWORTH, rp, rs )
sc_b, sc_a = scipy.signal.butter(order, [freq, freqStop], btype='bandpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'BandPass Butterworth order %d' % order)


### Bessel
# Test with even order
order = 4

rc = ricaudio.BandPass( order, freq, freqStop, ricaudio.BESSEL, rp, rs )
sc_b, sc_a = scipy.signal.bessel(order, [freq, freqStop], btype='bandpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'BandPass Bessel order %d' % order)


# Test with odd order
order = 5

rc = ricaudio.BandPass( order, freq, freqStop, ricaudio.BESSEL, rp, rs )
sc_b, sc_a = scipy.signal.bessel(order, [freq, freqStop], btype='bandpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'BandPass Bessel order %d' % order)


if plot:
    import pylab
    pylab.show()
