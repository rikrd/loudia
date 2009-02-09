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

a = [2, 0]

rr = ricaudio.poly( a )
rs = scipy.poly( a[:-1] )

print 'poly, plus zero:', scipy.allclose(rr[0,:-1], rs)


# Test the zpk <--> coeffs functions
z = scipy.array([1, -0.96, 0.80], dtype = 'f4')
p = scipy.array([1, 0.5, 0.5], dtype = 'f4')
k = 1.2

ra, rb = ricaudio.zpkToCoeffs(z, p, k)
sa, sb = scipy.signal.zpk2tf(z, p, k)

print 'zpkToCoeffs: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)


# Test the lpTolp functions
a = scipy.array([1, -0.96, 0.80], dtype = 'f4')
b = scipy.array([1, 0.5, 0.5], dtype = 'f4')

freq = 0.23

rb, ra = ricaudio.lowPassToLowPass(b, a, freq)
sb, sa = scipy.signal.lp2lp(b, a, freq)

print 'lpTolp: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)

# Test the comb function
rc = ricaudio.comb(5, 3)
sc = scipy.comb(5, 3)
print 'comb: ', round(sc) == round(rc)

# Test the bilinear function
a = scipy.array([10, -0.96, 0.80], dtype = 'f4')
b = scipy.array([156, 0.5, 0.5], dtype = 'f4')

rb, ra = ricaudio.bilinear(b, a, 1.0)
sb, sa = scipy.signal.bilinear(b, a)

# The ricaudio bilinear function does not return the coefficients normalized
rb = rb / ra[:,0]
ra = ra / ra[:,0]

print 'bilinear: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)
