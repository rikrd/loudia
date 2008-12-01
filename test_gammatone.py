#!/usr/bin/env python

# Create input
import scipy
import scipy.signal

samplerate = 44100
size = 4096*1
a1 = scipy.array(scipy.random.random((size, 1)), dtype='f4')
a2 = scipy.array(scipy.random.random((size, 1)), dtype='f4')

# Setup the Gammatone parameters and coefficients # --------------------- #
numFilters = 30

# Frequency parameters
lowFreq = 100.0
highFreq = samplerate / 2.0

# Resolution parameters
c = 9.26449
d = 24.7

def freqs_b_by_freqs(numFilters, lowFreq, highFreq, c, d, order = 1):
    # Find the center frequencies of the filters from the begin and end frequencies
    EarQ = c
    minBW = d
    vec = scipy.arange(numFilters, 0, -1)
    freqs = -(EarQ*minBW) + scipy.exp(vec*(-scipy.log(highFreq + EarQ*minBW) + scipy.log(lowFreq + EarQ*minBW))/numFilters) * (highFreq + EarQ*minBW);
    
    ERB = scipy.power((scipy.power((freqs/EarQ), order) + scipy.power(minBW, order)), (1.0 / order))
    B = 1.019 * 2.0 * scipy.pi * ERB
    
    return (freqs, B)
    
def filter_coeffs_gammatone(freqs, B, samplerate):
    T = 1 / float(samplerate)
    A0 = T
    A2 = 0.0
    B0 = 1.0
    B1 = -2.0 * scipy.cos(2.0 * freqs * scipy.pi * T) / scipy.exp(B * T)
    B2 = scipy.exp(-2.0 * B * T)
    
    A11 = -(2.0 * T * scipy.cos(2.0 * freqs * scipy.pi * T) / scipy.exp(B * T) + 2.0 * \
            scipy.sqrt(3+2.0**1.5) * T * scipy.sin(2.0 * freqs * scipy.pi * T) /  scipy.exp(B * T))/2.0
    
    A12 = -(2.0 * T * scipy.cos(2.0 * freqs * scipy.pi * T) / scipy.exp(B * T) - 2.0 * \
            scipy.sqrt(3+2.0**1.5) * T * scipy.sin(2.0 * freqs * scipy.pi * T) /  scipy.exp(B * T))/2.0
    
    A13 = -(2.0 * T * scipy.cos(2.0 * freqs * scipy.pi * T) / scipy.exp(B * T) + 2.0 * \
            scipy.sqrt(3-2.0**1.5) * T * scipy.sin(2.0 * freqs * scipy.pi * T) /  scipy.exp(B * T))/2.0
    
    A14 = -(2.0 * T * scipy.cos(2.0 * freqs * scipy.pi * T) / scipy.exp(B * T) - 2.0 * \
            scipy.sqrt(3-2.0**1.5) * T * scipy.sin(2.0 * freqs * scipy.pi * T) /  scipy.exp(B * T))/2.0
    
    gain = abs((-2.0 * scipy.exp(4 * 1j * freqs * scipy.pi * T) * T + \
                2.0 * scipy.exp(-(B * T) + 2.0 * 1j * freqs * scipy.pi * T) * T * \
                (scipy.cos(2.0 * freqs * scipy.pi * T) - scipy.sqrt(3 - 2.0**(3/2.0)) *  \
                 scipy.sin(2.0 * freqs * scipy.pi * T)))  *  \
               (-2.0 * scipy.exp(4 * 1j * freqs * scipy.pi * T) * T + \
                2.0 * scipy.exp(-(B * T) + 2.0 * 1j * freqs * scipy.pi * T) * T *  \
                (scipy.cos(2.0 * freqs * scipy.pi * T) + scipy.sqrt(3 - 2.0**(3/2.0))  *  \
                 scipy.sin(2.0 * freqs * scipy.pi * T))) *  \
               (-2.0 * scipy.exp(4 * 1j * freqs * scipy.pi * T) * T + \
                2.0 * scipy.exp(-(B * T) + 2.0 * 1j * freqs * scipy.pi * T) * T *  \
                (scipy.cos(2.0 * freqs * scipy.pi * T) - \
                 scipy.sqrt(3 + 2.0**(3/2.0)) * scipy.sin(2.0 * freqs * scipy.pi * T)))  *  \
               (-2.0 * scipy.exp(4 * 1j * freqs * scipy.pi * T) * T + 2.0 * scipy.exp(-(B * T) + 2.0 * 1j * freqs * scipy.pi * T) * T *  \
                (scipy.cos(2.0 * freqs * scipy.pi * T) + scipy.sqrt(3 + 2.0**(3/2.0)) * scipy.sin(2.0 * freqs * scipy.pi * T)))  /  \
               (-2.0  /  scipy.exp(2.0 * B * T) - 2.0 * scipy.exp(4 * 1j * freqs * scipy.pi * T) +  \
                2.0 * (1 + scipy.exp(4 * 1j * freqs * scipy.pi * T)) / scipy.exp(B * T))**4)
    
    allfilts = scipy.ones(len(freqs))
    fcoefs = (A0 * allfilts, A11, A12, A13, A14, A2 * allfilts, B0 * allfilts, B1, B2, gain)
    return fcoefs
        
freqs, B = freqs_b_by_freqs(numFilters, lowFreq, highFreq, c, d)
(A0, A11, A12, A13, A14, A2, B0, B1, B2, gain) = filter_coeffs_gammatone(freqs, B, samplerate)
# -------------------------------------------------------- #


# CRicaudio's solution # --------------------------------- #
import pycricaudio

coeffsA1 = scipy.array(scipy.vstack((A0, A11, A2)) / gain, dtype = 'f4')
coeffsA2 = scipy.array(scipy.vstack((A0, A12, A2)), dtype = 'f4')
coeffsA3 = scipy.array(scipy.vstack((A0, A13, A2)), dtype = 'f4')
coeffsA4 = scipy.array(scipy.vstack((A0, A14, A2)), dtype = 'f4')

coeffsB = scipy.array(scipy.vstack((B0, B1, B2)), dtype = 'f4')

f1 = pycricaudio.Filter(coeffsB, coeffsA1, samplerate, numFilters)
f2 = pycricaudio.Filter(coeffsB, coeffsA2, samplerate, numFilters)
f3 = pycricaudio.Filter(coeffsB, coeffsA3, samplerate, numFilters)
f4 = pycricaudio.Filter(coeffsB, coeffsA4, samplerate, numFilters)

b1 = f4.process(f3.process(f2.process(f1.process(a1))))
b2 = f4.process(f3.process(f2.process(f1.process(a2))))

print 'CRicaudio: Done'
# -------------------------------------------------------- #

# Ricaudio's solution # ---------------------------------- #    
zi = scipy.zeros((0,0))

def filterbank_compute(samples):
    v = samples
    x = scipy.resize(v, (gain.shape[0], v.shape[0]))
    
    if zi.shape[0] != gain.shape[0]:
        zi.resize((max(gain.shape[0], gain.shape[0]), 4 , 2))
        
    def filt(x):
        coeffsA1 = scipy.array([A0[row[0]] / gain[row[0]],
                                A11[row[0]]/ gain[row[0]],
                                A2[row[0]] / gain[row[0]]], dtype = 'f4')

        b = scipy.array([B0[row[0]], B1[row[0]], B2[row[0]]])

        y1, zi[row[0],0,:] = scipy.signal.lfilter(coeffsA1,
                                                  b,
                                                  x, zi = zi[row[0],0,:])
        
        y2, zi[row[0],1,:] = scipy.signal.lfilter([A0[row[0]],
                                                   A12[row[0]],
                                                   A2[row[0]]],
                                                  b,
                                                  y1, zi = zi[row[0],1,:])
        
        y3, zi[row[0],2,:] = scipy.signal.lfilter([A0[row[0]],
                                                   A13[row[0]],
                                                   A2[row[0]]],
                                                  b,
                                                  y2, zi = zi[row[0],2,:])
        
        y4, zi[row[0],3,:] = scipy.signal.lfilter([A0[row[0]],
                                                   A14[row[0]],
                                                   A2[row[0]]],
                                                  b,
                                                  y3, zi = zi[row[0],3,:])
        row[0] += 1
        return y4

    row = [0]
    y = scipy.apply_along_axis(filt, 1, x)
    return y.T 

c1 = filterbank_compute(a1)
c2 = filterbank_compute(a2)

print 'Ricaudio: Done'
# -------------------------------------------------------- #

print 'Difference:'
print '  %f' % ((abs(b1-c1)).sum() / abs(c1).sum(),)
print '  %f' % ((abs(b2-c2)).sum() / abs(c2).sum(),)

print scipy.allclose(c1, b1, rtol = 1e2)
print scipy.allclose(c2, b2, rtol = 1e2)
