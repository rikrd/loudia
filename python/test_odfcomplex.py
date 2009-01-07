#!/usr/bin/env python

import ricaudio
import scipy

a = [[1  ,  2,  6,  4],
     [1.2,  6,  2,  4],
     [1  ,  2,  6,  4],
     [30  , 20, 60, 40]]

d = ricaudio.ODFComplex(4, 16)
b1 = d.process(a)

print b1
