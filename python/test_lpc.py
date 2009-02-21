#!/usr/bin/env python

import ricaudio
import scipy

d = ricaudio.LPC(9, 9)
a1 = scipy.array([[1,2,3,4,5,6,7,8,9]], dtype = 'f4')
a2 = scipy.array([[0,0,0,0,0,0,0,0,0]], dtype = 'f4')

coeffs1, reflection1, error1 = d.process(a1)
coeffs2, reflection2, error2 = d.process(a2)

exp_coeffs1 = scipy.array([[ 1.        , -0.89023799,  0.01518036,  0.01562015,  0.01640711,
                             0.01755869,  0.01910042,  0.02106662, -0.00439638]], dtype=scipy.float32)

exp_reflection1 = scipy.array([[-0.84210527,  0.07365081,  0.06710767,  0.05878957,  0.04811161,
                                0.03444673,  0.01715312, -0.00439638]], dtype=scipy.float32)

exp_error1 = scipy.array([[ 81.4784317]], dtype=scipy.float32)



exp_coeffs2 = scipy.array([[1,0,0,0,0,0,0,0,0]], dtype=scipy.float32)

exp_reflection2 = scipy.array([[0,0,0,0,0,0,0,0]], dtype=scipy.float32)

exp_error2 = scipy.array([[0]], dtype=scipy.float32)


print scipy.allclose(coeffs1, exp_coeffs1)
print scipy.allclose(reflection1, exp_reflection1)
print scipy.allclose(error1, exp_error1)

print scipy.allclose(coeffs2, exp_coeffs2)
print scipy.allclose(reflection2, exp_reflection2)
print scipy.allclose(error2, exp_error2)
