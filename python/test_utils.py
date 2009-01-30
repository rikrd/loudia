#!/usr/bin/env python

import ricaudio
import scipy

# Test the magToDb <-> dbToMag
i = scipy.array([[1,2,3,4,5,6]], dtype = 'f4')
a = ricaudio.magToDb(i)
o = ricaudio.dbToMag(a)

print scipy.allclose(i, o)

# Test the transform of windows
transf = ricaudio.hammingTransform(24, 10, 1024, 4096)

print max(transf)
