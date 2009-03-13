#!/usr/bin/env python

import ricaudio
import scipy
import scipy.signal
from common import *

plot = True

atol = 1e-5

freq = 0.3
fs = 8000

### Chebyshev I
# Test with even order
rp = 0.05
order = 4

rc = ricaudio.HighPass( order, freq, rp, ricaudio.CHEBYSHEVI )
sc_b, sc_a = scipy.signal.cheby1(order, rp, freq, btype='highpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'HighPass Chebyshev I order %d' % order)

# Test with odd order
order = 5

rc = ricaudio.HighPass( order, freq, rp, ricaudio.CHEBYSHEVI )
sc_b, sc_a = scipy.signal.cheby1( order, rp, freq, btype='highpass', analog=0, output='ba' )

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'HighPass Chebyshev I order %d' % order)


### Chebyshev II
# Test with even order
rp = 40
order = 4

rc = ricaudio.HighPass( order, freq, rp, ricaudio.CHEBYSHEVII )
sc_b, sc_a = scipy.signal.cheby2(order, rp, freq, btype='highpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'HighPass Chebyshev II order %d' % order)


# Test with odd order
order = 5

rc = ricaudio.HighPass( order, freq, rp, ricaudio.CHEBYSHEVII )
sc_b, sc_a = scipy.signal.cheby2(order, rp, freq, btype='highpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'HighPass Chebyshev II order %d' % order)


### Butterworth
# Test with even order
rp = 40
order = 4

rc = ricaudio.HighPass( order, freq, rp, ricaudio.BUTTERWORTH )
sc_b, sc_a = scipy.signal.butter(order, freq, btype='highpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'HighPass Butterworth order %d' % order)


# Test with odd order
order = 5

rc = ricaudio.HighPass( order, freq, rp, ricaudio.BUTTERWORTH )
sc_b, sc_a = scipy.signal.butter(order, freq, btype='highpass', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    plotFreqz(rc_b, rc_a, title = 'HighPass Butterworth order %d' % order)

if plot:
    import pylab
    pylab.show()
