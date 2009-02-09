#!/usr/bin/env python

# Create input
import scipy
import scipy.signal
import ricaudio

order = 8
freq = 0.3
rp = 0.05
fs = 8000

rc = ricaudio.Chebyshev( 1, order, freq, rp, 8000 )
sc_b, sc_a = scipy.signal.cheby1(order, rp, freq, btype='low', analog=0, output='ba')

print scipy.allclose(sc_b, rc.b().T) and scipy.allclose(sc_a, rc.a().T)


