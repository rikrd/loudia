#!/usr/bin/env python

import ricaudio
import scipy

d = ricaudio.LPC(9, 9)
a = scipy.array([[1,2,3,4,5,6,7,8,9]], dtype = 'f4')
b = d.process(a)

print b
