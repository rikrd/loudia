#!/usr/bin/env python

import ricaudio
import scipy
import scipy.signal

# Test the magToDb <-> dbToMag
i = scipy.array([[1,2,3,4,5,6]], dtype = 'f4')
a = ricaudio.magToDb(i)
o = ricaudio.dbToMag(a)

print 'magToDb:', scipy.allclose(i, o)

# Test the transform of windows
transf = ricaudio.hammingTransform(24, 10, 1024, 4096)


# Test the poly function
a = [1, 2, 3, 4, 5]

rr = ricaudio.poly( a )
rs = scipy.poly( a )

print 'poly:', scipy.allclose(rr, rs)


# Test the zpk <--> coeffs functions
z = scipy.array([1, -0.96, 0.80], dtype = 'f4')
p = scipy.array([1, 0.5, 0.5], dtype = 'f4')
k = 1.2

ra, rb = ricaudio.zpkToCoeffs(z, p, k)
sa, sb = scipy.signal.zpk2tf(z, p, k)

print 'zpkToCoeffs: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)
