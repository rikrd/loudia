#!/usr/bin/env python

# Create input
import scipy
import loudia

a = scipy.array([1, 2, 3, 4, 5])

rr = loudia.poly( a )
rs = scipy.poly( a )

print rr[0]
print rs
print scipy.allclose(rr[0], rs, rtol = 1e1)


