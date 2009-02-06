#!/usr/bin/env python

# Create input
import scipy
import ricaudio

plot = True
frameSize = 256
fftSize = 512
samplerate = 8000

a_zeros = scipy.array(scipy.zeros((1, frameSize)), dtype='f4')
a_ones = scipy.array(scipy.ones((1, frameSize)), dtype='f4')
a_random = scipy.array(scipy.random.random((1, frameSize)), dtype='f4')
a_sine = scipy.array(scipy.cos(2 * scipy.pi * 440 * scipy.arange(frameSize) / samplerate + scipy.pi/4.0), dtype='f4')
a_sine = a_sine.reshape((1, a_sine.shape[0]))

# Ricaudio's solution # --------------------------------- #
w = ricaudio.Window(frameSize, ricaudio.Window.HAMMING)
m = ricaudio.FFT(frameSize, fftSize, True)
n = ricaudio.IFFT(fftSize, frameSize, True)

r_zeros = n.process(m.process(w.process(a_zeros)))
r_ones = n.process(m.process(w.process(a_ones)))
r_random = n.process(m.process(w.process(a_random)))
r_sine = n.process(m.process(w.process(a_sine)))
# -------------------------------------------------------- #

x_zeros = w.process(a_zeros)
x_ones = w.process(a_ones)
x_random = w.process(a_random)
x_sine = w.process(a_sine)

atol = 1e-6

print scipy.allclose(r_zeros, x_zeros, atol = atol)
print scipy.allclose(r_ones, x_ones, atol = atol)
print scipy.allclose(r_random, x_random, atol = atol)
print scipy.allclose(r_sine, x_sine, atol = atol)

if plot:
    import pylab
    pylab.subplot(211)
    pylab.hold(True)
    pylab.plot(r_sine.T, label = 'Ricaudio')
    pylab.plot(x_sine.T, label = 'Expected')

    pylab.legend()

    
    pylab.subplot(212)
    pylab.hold(True)
    pylab.plot(r_sine.T - x_sine.T, label = 'Difference')
    
    pylab.legend()
    pylab.show()
