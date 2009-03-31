#!/usr/bin/env python

import loudia
import scipy
import scipy.signal

# Test the magToDb <-> dbToMag
i = scipy.arange(1, 7, dtype = 'f4')
a = loudia.magToDb(i)
o = loudia.dbToMag(a)

print 'magToDb/dbToMag:', scipy.allclose(i, o)

# Test the transform of windows
transf = loudia.hammingTransform(24, 10, 1024, 4096)


# Test the poly function
a = scipy.array([1, 2, 3, 4, 5])

rr = loudia.poly( a )
rs = scipy.poly( a )

print 'poly:', scipy.allclose(rr, rs)

a = scipy.array([2, 0])

rr = loudia.poly( a )
rs = scipy.poly( a[:-1] )

print 'poly, plus zero:', scipy.allclose(rr[0,:-1], rs)


# Test the zpk <--> coeffs functions
z = scipy.array([1, -0.96, 0.80])
p = scipy.array([1, 0.5, 0.5])
k = 1.2

ra, rb = loudia.zpkToCoeffs(z, p, k)
sa, sb = scipy.signal.zpk2tf(z, p, k)

print 'zpkToCoeffs: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)


# Test the lpTolp functions
a = scipy.array([1, -0.96, 0.80])
b = scipy.array([1, 0.5, 0.5])

freq = 0.23

rb, ra = loudia.lowPassToLowPass(b, a, freq)
sb, sa = scipy.signal.lp2lp(b, a, freq)

print 'lpTolp: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)

# Test the comb function
rc = loudia.combination(5, 3)
sc = scipy.comb(5, 3)
print 'combination: ', round(sc) == round(rc)

# Test the bilinear function
a = scipy.array([10, -0.96, 0.80])
b = scipy.array([156, 0.5, 0.5])

rb, ra = loudia.bilinear(b, a, 1.0)
sb, sa = scipy.signal.bilinear(b, a)

# The loudia bilinear function does not return the coefficients normalized
rb = rb / ra[:,0]
ra = ra / ra[:,0]

print 'bilinear: ', scipy.allclose(rb, sb) and scipy.allclose(ra, sa)

# Test the nextPowerOf2 function
print 'nextPowerOf2(0): ', loudia.nextPowerOf2(0) == 2
print 'nextPowerOf2(0, 3): ', loudia.nextPowerOf2(0, 3) == 16
print 'nextPowerOf2(94): ', loudia.nextPowerOf2(94) == 128
print 'nextPowerOf2(94, 1): ', loudia.nextPowerOf2(94, 1) == 256
