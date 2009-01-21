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
d5 = ricaudio.ODF(8, ricaudio.ODF.NORM_WEIGHTED_PHASE_DEVIATION)
d6 = ricaudio.ODF(8, ricaudio.ODF.SPECTRAL_FLUX)
d7 = ricaudio.ODF(8, ricaudio.ODF.MODIFIED_KULLBACK_LIEBLER)

b1 = d1.process(a)
b2 = d2.process(a)
b3 = d3.process(a)
b4 = d4.process(a)
b5 = d5.process(a)
b6 = d6.process(a)
b7 = d7.process(a)

print b1[:,0], ': complex domain'
print b2[:,0], ': rectified complex domain'
print b3[:,0], ': phase deviation'
print b4[:,0], ': weighted phase deviation'
print b5[:,0], ': normalized weighted phase deviation'
print b6[:,0], ': spectral flux'
print b7[:,0], ': modified kullback-liebler'
