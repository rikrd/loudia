#!/usr/bin/env python

import ricaudio
import scipy

a = [[1  , 2, 6, 4],
     [1.2, 6, 2, 4]]

d = ricaudio.Unwrap(4)
b1 = d.process(a)

# Use a different unwarp since scipy's is 2x slower
def fastunwrap(thetaArray, discont = scipy.pi):
    # takes an array of theta values
    # returns the data in unwrapped form (unwrapping over the axis == 1)
    diff = scipy.zeros_like(thetaArray)
    diff[1:,:] = scipy.diff(thetaArray, axis = 0)
    upSteps = diff > discont
    downSteps = diff < -discont
    shift = scipy.cumsum(upSteps, axis = 0) - scipy.cumsum(downSteps, axis = 0)
    return thetaArray - 2.0*discont*shift

b2 = fastunwrap(a)

print b1
print b2

print scipy.allclose(b1, b2)
