#!/usr/bin/env python

import loudia
import scipy

d1 = loudia.LPC(9, 9)
d2 = loudia.LPCResidual(9)
a = scipy.arange(1,10)
coeffs, reflection, error = d1.process(a)
residual = d2.process(a, coeffs)

exp_coeffs = [ 1.        , -0.89023799,  0.01518036,  0.01562015,  0.01640711,
               0.01755869,  0.01910042,  0.02106662, -0.00439638]

exp_reflection = [-0.84210527,  0.07365081,  0.06710767,  0.05878957,  0.04811161,
                  0.03444673,  0.01715312, -0.00439638]

exp_error = [ 81.4784317]


print scipy.allclose(coeffs, exp_coeffs)
print scipy.allclose(reflection, exp_reflection)
print scipy.allclose(error, exp_error)
