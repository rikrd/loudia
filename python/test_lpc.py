#!/usr/bin/env python

import ricaudio

d = ricaudio.LPC(9, 9)
b = d.process([1,2,3,4,5,6,7,8,9])

print b
