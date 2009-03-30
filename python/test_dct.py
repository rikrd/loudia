#!/usr/bin/env python

import loudia
import scipy

d = loudia.DCT(9, 9, True, loudia.DCT.OCTAVE)
b = d.process( scipy.array( [1,2,3,4,5,6,7,8,9] ) )

print b
