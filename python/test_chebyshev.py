#!/usr/bin/env python

# Create input
import scipy
import scipy.signal
import ricaudio

plot = True

atol = 1e-5

freq = 0.3
fs = 8000

# Test type I with even order
rp = 0.05
order = 4

rc = ricaudio.Chebyshev( order, freq, rp, 1 )
sc_b, sc_a = scipy.signal.cheby1(order, rp, freq, btype='low', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    # Create the omega array
    npoints = 1000
    w = scipy.arange(-scipy.pi, scipy.pi, 2*scipy.pi/(npoints), dtype = 'f4')

    # Calculate the frequency response
    d = ricaudio.freqz(rc_b.T, rc_a.T, w)

    import pylab

    pylab.subplot(2,1,1)
    pylab.plot(w, abs(d[:,0]))
    pylab.title('Type I with order %d \n Magnitude of the Frequency Response' % order)
    
    pylab.subplot(2,1,2)
    pylab.plot(w, scipy.angle(d[:,0]))
    pylab.title('Angle of the Frequency Response')
    
    pylab.show()


# Test with type I with odd order
order = 11

rc = ricaudio.Chebyshev( order, freq, rp, 1, ricaudio.Chebyshev.I )
sc_b, sc_a = scipy.signal.cheby1( order, rp, freq, btype='low', analog=0, output='ba' )

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    # Create the omega array
    npoints = 1000
    w = scipy.arange(-scipy.pi, scipy.pi, 2*scipy.pi/(npoints), dtype = 'f4')

    # Calculate the frequency response
    d = ricaudio.freqz(rc_b.T, rc_a.T, w)

    import pylab

    pylab.subplot(2,1,1)
    pylab.plot(w, abs(d[:,0]))
    pylab.title('Type I with order %d \n Magnitude of the Frequency Response' % order)
    
    pylab.subplot(2,1,2)
    pylab.plot(w, scipy.angle(d[:,0]))
    pylab.title('Angle of the Frequency Response')
    
    pylab.show()


# Test type II with even order
rc = ricaudio.Chebyshev( order, freq, rp, 1, ricaudio.Chebyshev.II )
sc_b, sc_a = scipy.signal.cheby2(order, rp, freq, btype='low', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    # Create the omega array
    npoints = 1000
    w = scipy.arange(-scipy.pi, scipy.pi, 2*scipy.pi/(npoints), dtype = 'f4')

    # Calculate the frequency response
    d = ricaudio.freqz(rc_b.T, rc_a.T, w)

    import pylab

    pylab.subplot(2,1,1)
    pylab.plot(w, abs(d[:,0]))
    pylab.title('Type II with order %d \n Magnitude of the Frequency Response' % order)
    
    pylab.subplot(2,1,2)
    pylab.plot(w, scipy.angle(d[:,0]))
    pylab.title('Angle of the Frequency Response')
    
    pylab.show()


# Test with type II with odd order
rp = 5
order = 11

rc = ricaudio.Chebyshev( order, freq, rp, 1, ricaudio.Chebyshev.II )
sc_b, sc_a = scipy.signal.cheby2(order, rp, freq, btype='low', analog=0, output='ba')

rc_b = rc.b().T
rc_a = rc.a().T

print scipy.allclose(sc_b, rc_b, atol = atol) and scipy.allclose(sc_a, rc_a, atol = atol)

if plot:
    # Create the omega array
    npoints = 1000
    w = scipy.arange(-scipy.pi, scipy.pi, 2*scipy.pi/(npoints), dtype = 'f4')

    # Calculate the frequency response
    d = ricaudio.freqz(rc_b.T, rc_a.T, w)

    import pylab

    pylab.subplot(2,1,1)
    pylab.plot(w, abs(d[:,0]))
    pylab.title('Type II with order %d \n Magnitude of the Frequency Response' % order)
    
    pylab.subplot(2,1,2)
    pylab.plot(w, scipy.angle(d[:,0]))
    pylab.title('Angle of the Frequency Response')
    
    pylab.show()
