#!/usr/bin/env python

# Create input
import scipy
import scipy.signal
import ricaudio

order = 12
freq = 0.3
rp = -60
fs = 8000

rc = ricaudio.Chebyshev( 1, order, freq, rp, 8000 )
sc = cheby1(order, rp, freq, btype='low', analog=0, output='ba')

print sc
print rc.a()
print rc.b()

#print scipy.allclose(rr[0], rs, rtol = 1e1)


