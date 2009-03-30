#!/usr/bin/env python

import loudia
import scipy

a = [[10,  20,  6,  4],
     [100,  10,  13,  10],
     [10,  100,  20,  0],
     [200,  10,  13,  10],
     [100,  210,  313,  510]]

d1 = loudia.ODF(8, loudia.ODF.COMPLEX_DOMAIN)
d2 = loudia.ODF(8, loudia.ODF.RECTIFIED_COMPLEX_DOMAIN)
d3 = loudia.ODF(8, loudia.ODF.PHASE_DEVIATION)
d4 = loudia.ODF(8, loudia.ODF.WEIGHTED_PHASE_DEVIATION)
d5 = loudia.ODF(8, loudia.ODF.NORM_WEIGHTED_PHASE_DEVIATION)
d6 = loudia.ODF(8, loudia.ODF.SPECTRAL_FLUX)
d7 = loudia.ODF(8, loudia.ODF.MODIFIED_KULLBACK_LIEBLER)
d8 = loudia.ODF(8, loudia.ODF.HIGH_FREQUENCY_CONTENT)

b1 = d1.process(a)
b2 = d2.process(a)
b3 = d3.process(a)
b4 = d4.process(a)
b5 = d5.process(a)
b6 = d6.process(a)
b7 = d7.process(a)
b8 = d8.process(a)

print b1[:,0], ': complex domain'
print b2[:,0], ': rectified complex domain'
print b3[:,0], ': phase deviation'
print b4[:,0], ': weighted phase deviation'
print b5[:,0], ': normalized weighted phase deviation'
print b6[:,0], ': spectral flux'
print b7[:,0], ': modified kullback-liebler'
print b8[:,0], ': high frequency content'
