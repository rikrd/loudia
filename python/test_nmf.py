#!/usr/bin/env python

import ricaudio
import scipy
import pylab

plot = False

d = ricaudio.NMF(8, 2, 100, 1e-7)

a = scipy.zeros((14, 8))
a[:4, 2] = 1
a[5:9, 6] = 1
a[10:, 2] = 1

expectedComponents = [[ 0., 0., 0.        , 0., 0., 0., 1.82457483, 0. ],
                      [ 0., 0., 2.91204405, 0., 0., 0., 0.        , 0. ]]

expectedGains = [[ 0. , 0.35355338],
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
                 [ 0. , 0.35355338]]

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
