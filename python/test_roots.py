#!/usr/bin/env python

# Create input
import scipy
import ricaudio

channels = 1
order = 10
rippleDB = 1
samplerate = 44100

poly = [1, 2, 3, 4, 5]

c = ricaudio.Chebyshev(channels, order, rippleDB, samplerate)
rr = c.roots( poly )

rs = scipy.roots( poly )


print rr[0]
print rs
print scipy.allclose(rr[0], rs, rtol = 1e1)


