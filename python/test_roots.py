#!/usr/bin/env python

# Create input
import scipy
import ricaudio

poly = [1, 2, 3, 4, 5]

rr = ricaudio.roots( poly )

rs = scipy.flipud(scipy.roots( poly ))


print rr[:,0]
print rs
print scipy.allclose(rr[0], rs, rtol = 1e1)


