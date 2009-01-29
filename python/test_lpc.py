#!/usr/bin/env python

import ricaudio
import scipy

d = ricaudio.LPC(9, 9)
a = scipy.array([[1,2,3,4,5,6,7,8,9]], dtype = 'f4')
coeffs, reflection, error = d.process(a)


exp_coeffs = scipy.array([[ 1.        , -0.89023799,  0.01518036,  0.01562015,  0.01640711,
                            0.01755869,  0.01910042,  0.02106662, -0.00439638]], dtype=scipy.float32)

exp_reflection = scipy.array([[-0.84210527,  0.07365081,  0.06710767,  0.05878957,  0.04811161,
                               0.03444673,  0.01715312, -0.00439638]], dtype=scipy.float32)

exp_error = scipy.array([[ 81.4784317]], dtype=scipy.float32)


print scipy.allclose(coeffs, exp_coeffs)
print scipy.allclose(reflection, exp_reflection)
print scipy.allclose(error, exp_error)
