#!/usr/bin/env python

# Create input
import scipy
import loudia

a = scipy.array([1, 2, 3, 4, 5])
b = scipy.array([2, 4, 5])

rr = loudia.convolve( a, b )
rs = scipy.convolve( a, b )

print rr[0]
print rs
print scipy.allclose(rr[0], rs)


