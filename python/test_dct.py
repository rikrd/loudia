#!/usr/bin/env python

import ricaudio
import scipy

d = ricaudio.DCT(9, 9, True, ricaudio.DCT.OCTAVE)
b = d.process( scipy.array( [1,2,3,4,5,6,7,8,9] ) )

print b
