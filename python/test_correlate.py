#!/usr/bin/env python

import ricaudio
import scipy

a = scipy.array([[1,2,3,4,5,6,7,8,9]], dtype = 'f4')
b = scipy.array([[1,2,3,4,5,6,7,8,9]], dtype = 'f4')

r = ricaudio.correlate(a, b)
s = scipy.correlate(a[0,:], b[0,:], 'full')
print r
print scipy.allclose(r[0,:], s)

r = ricaudio.correlate(b, a)
s = scipy.correlate(b[0,:], a[0,:], 'full')

print scipy.allclose(r[0,:], s)

