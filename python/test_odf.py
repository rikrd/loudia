#!/usr/bin/env python

import ricaudio
import scipy

a = [[10,  20,  6,  4],
     [100,  10,  13,  10],
     [10,  100,  20,  0]]

d1 = ricaudio.ODF(8, ricaudio.ODF.COMPLEX_DOMAIN)

d2 = ricaudio.ODF(8, ricaudio.ODF.RECTIFIED_COMPLEX_DOMAIN)

b1 = d1.process(a)
b2 = d2.process(a)

print b1[:,0]
print b2[:,0]
