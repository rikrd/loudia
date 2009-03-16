#!/usr/bin/env python

import ricaudio
import scipy
import scipy.signal

# Test the magToDb <-> dbToMag
i = scipy.arange(1, 7)
a = ricaudio.magToDb(i)
o = ricaudio.dbToMag(a)

print 'magToDb/dbToMag:', scipy.allclose(i, o)

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
z = [1, -0.96, 0.80]
p = [1, 0.5, 0.5]
k = 1.2

ra, rb = ricaudio.zpkToCoeffs(z, p, k)
sa, sb = scipy.signal.zpk2tf(z, p, k)

print 'zpkToCoeffs: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)


# Test the lpTolp functions
a = [1, -0.96, 0.80]
b = [1, 0.5, 0.5]

freq = 0.23

rb, ra = ricaudio.lowPassToLowPass(b, a, freq)
sb, sa = scipy.signal.lp2lp(b, a, freq)

print 'lpTolp: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)

# Test the comb function
rc = ricaudio.comb(5, 3)
sc = scipy.comb(5, 3)
print 'comb: ', round(sc) == round(rc)

# Test the bilinear function
a = [10, -0.96, 0.80]
b = [156, 0.5, 0.5]

rb, ra = ricaudio.bilinear(b, a, 1.0)
sb, sa = scipy.signal.bilinear(b, a)

# The ricaudio bilinear function does not return the coefficients normalized
rb = rb / ra[:,0]
ra = ra / ra[:,0]

print 'bilinear: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)

# Test the nextPowerOf2 function
print 'nextPowerOf2(0): ', ricaudio.nextPowerOf2(0) == 2
print 'nextPowerOf2(0, 3): ', ricaudio.nextPowerOf2(0, 3) == 16
print 'nextPowerOf2(94): ', ricaudio.nextPowerOf2(94) == 128
print 'nextPowerOf2(94, 1): ', ricaudio.nextPowerOf2(94, 1) == 256
