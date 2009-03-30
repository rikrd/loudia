#!/usr/bin/env python

import loudia
import scipy

rtol = 1e-5
atol = 1e-7

d = loudia.LPC(9, 9, 0.99)
a1 = scipy.arange(1, 10)
a2 = scipy.zeros(9)

coeffs1, reflection1, error1 = d.process(a1)
coeffs2, reflection2, error2 = d.process(a2)


d = loudia.LPC(9, 9)
a1 = scipy.arange(1,10)
a2 = scipy.zeros(9)

coeffs1, reflection1, error1 = d.process(a1)
coeffs2, reflection2, error2 = d.process(a2)



exp_coeffs1 = [ 1.        , -0.89023799,  0.01518036,  0.01562015,  0.01640711,
                0.01755869,  0.01910042,  0.02106662, -0.00439638]

exp_reflection1 = [-0.84210527,  0.07365081,  0.06710767,  0.05878957,  0.04811161,
                   0.03444673,  0.01715312, -0.00439638]

exp_error1 = [ 81.4784317]



exp_coeffs2 = [1,0,0,0,0,0,0,0,0]
exp_reflection2 = scipy.zeros(8)
exp_error2 = scipy.zeros(1)


print scipy.allclose(coeffs1, exp_coeffs1, rtol = rtol, atol = atol)
print scipy.allclose(reflection1, exp_reflection1, rtol = rtol, atol = atol)
print scipy.allclose(error1, exp_error1, rtol = rtol, atol = atol)

print scipy.allclose(coeffs2, exp_coeffs2, rtol = rtol, atol = atol)
print scipy.allclose(reflection2, exp_reflection2, rtol = rtol, atol = atol)
print scipy.allclose(error2, exp_error2, rtol = rtol, atol = atol)
