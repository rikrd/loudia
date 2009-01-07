#!/usr/bin/env python

import ricaudio

d = ricaudio.Unwrap(9)
b = d.process([[1,2,3,4,5,6,7,8,9],
               [1+3j,2,3+4j,4,5,6,7,8,9]])

print b
