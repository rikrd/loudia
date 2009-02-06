#!/usr/bin/env python

import ricaudio
import scipy
import pylab

plot = True

d = ricaudio.INMF(8, 2, 10, 0.5, 0.5, 100, 1, 1e-7)

a = scipy.zeros((14, 8), dtype = 'f4')
a[:4, 2] = 1
a[5:9, 6] = 1
a[10:, 2] = 1

expectedComponents = scipy.array([[ 0., 0., 0.        , 0., 0., 0., 1.82457483, 0. ],
                                  [ 0., 0., 2.91204405, 0., 0., 0., 0.        , 0. ]],
                                 dtype = 'f4')

expectedGains = scipy.array([[ 0. , 0.35355338],
                             [ 0. , 0.35355338],
                             [ 0. , 0.35355338],
                             [ 0. , 0.35355338],
                             [ 0. , 0.        ],
                             [ 0.5, 0.        ],
                             [ 0.5, 0.        ],
                             [ 0.5, 0.        ],
                             [ 0.5, 0.        ],
                             [ 0. , 0.        ],
                             [ 0. , 0.35355338],
                             [ 0. , 0.35355338],
                             [ 0. , 0.35355338],
                             [ 0. , 0.35355338]],
                            dtype = 'f4')

components, gains = d.process(a)

print 'components:', scipy.allclose(components, expectedComponents)
print 'gains:', scipy.allclose(gains, expectedGains)

if plot:
    pylab.figure()
    pylab.imshow(scipy.flipud(a.T))
    pylab.title("Input")
    
    pylab.figure()
    pylab.plot(components.T)
    pylab.title("Components")
    
    pylab.figure()
    pylab.plot(gains)
    pylab.title("Gains")
    
    pylab.show()
