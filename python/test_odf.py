#!/usr/bin/env python

import ricaudio
import scipy

a = [[10,  20,  6,  4],
     [100,  10,  13,  10],
     [10,  100,  20,  0]]

d1 = ricaudio.ODF(8, ricaudio.ODF.COMPLEX_DOMAIN)
d2 = ricaudio.ODF(8, ricaudio.ODF.RECTIFIED_COMPLEX_DOMAIN)
d3 = ricaudio.ODF(8, ricaudio.ODF.PHASE_DEVIATION)
d4 = ricaudio.ODF(8, ricaudio.ODF.WEIGHTED_PHASE_DEVIATION)
d5 = ricaudio.ODF(8, ricaudio.ODF.SPECTRAL_FLUX)

b1 = d1.process(a)
b2 = d2.process(a)

print b1[:,0]
print b2[:,0]
