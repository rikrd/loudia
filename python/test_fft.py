#!/usr/bin/env python

# Create input
import scipy
import loudia
import pylab

plot = False
frameSize = 121
fftSize = 512
sampleRate = 8000

a_zeros = scipy.zeros( frameSize )
a_ones = scipy.ones( frameSize )
a_random = scipy.random.random( frameSize )
a_sine = scipy.cos(2 * scipy.pi * 440 * scipy.arange( frameSize ) / sampleRate + scipy.pi/4.0)

# Loudia's solution # --------------------------------- #
w = loudia.Window(frameSize, loudia.Window.RECTANGULAR)
m = loudia.FFT(fftSize, False)

r_zeros = m.process(w.process(a_zeros))
r_ones = m.process(w.process(a_ones))
r_random = m.process(w.process(a_random))
r_sine = m.process(w.process(a_sine))
# -------------------------------------------------------- #


# Scipy solution # ---------------------------------- #
s_zeros = scipy.fft(a_zeros, fftSize)[:fftSize/2+1]
s_ones = scipy.fft(a_ones, fftSize)[:fftSize/2+1]
s_random = scipy.fft(a_random, fftSize)[:fftSize/2+1]
s_sine = scipy.fft(a_sine, fftSize)[:fftSize/2+1]
# -------------------------------------------------------- #

atol = 1e-5

print scipy.allclose(r_zeros, s_zeros, atol = atol)
print scipy.allclose(r_ones, s_ones, atol = atol)
print scipy.allclose(r_random, s_random, atol = atol)
print scipy.allclose(r_sine, s_sine, atol = atol)

if plot:
    r_abs = scipy.absolute(r_sine).T
    r_ang = scipy.angle(r_sine).T
    r_max = max(r_abs)
    
    s_abs = scipy.absolute(s_sine).T
    s_ang = scipy.angle(s_sine).T
    s_max = max(s_abs)

    pylab.subplot(211)
    pylab.hold(True)
    pylab.plot(r_abs, label = 'Loudia')
    pylab.plot(s_abs, label = 'Scipy')
    
    pylab.subplot(212)
    pylab.hold(True)
    pylab.plot(r_abs*r_ang/r_max, label = 'Loudia')
    pylab.plot(s_abs*s_ang/s_max, label = 'Scipy')
    
    pylab.legend()


fftSize = 1024

m = loudia.FFT(fftSize, False)
r_zeros = m.process(w.process(a_zeros))
r_ones = m.process(w.process(a_ones))
r_random = m.process(w.process(a_random))
r_sine = m.process(w.process(a_sine))

s_zeros = scipy.fft(a_zeros, fftSize)[:fftSize/2+1]
s_ones = scipy.fft(a_ones, fftSize)[:fftSize/2+1]
s_random = scipy.fft(a_random, fftSize)[:fftSize/2+1]
s_sine = scipy.fft(a_sine, fftSize)[:fftSize/2+1]

atol = 1e-5

print scipy.allclose(r_zeros, s_zeros, atol = atol)
print scipy.allclose(r_ones, s_ones, atol = atol)
print scipy.allclose(r_random, s_random, atol = atol)
print scipy.allclose(r_sine, s_sine, atol = atol)

if plot:
    pylab.show()
