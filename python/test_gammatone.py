#!/usr/bin/env python

# Create input
import scipy
import scipy.signal

# Plotting parameters
plotAngle = False
plotColor = False

# Setup the Gammatone parameters and coefficients # --------------------- #
numFilters = 30
sampleRate = 44100

# Frequency parameters
lowFreq = 100.0
highFreq = sampleRate / 2.0

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

def filter_coeffs_gammatone(freqs, B, sampleRate):
    T = 1 / float(sampleRate)
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
(B0, B11, B12, B13, B14, B2, A0, A1, A2, gain) = filter_coeffs_gammatone(freqs, B, sampleRate)
# -------------------------------------------------------- #


import loudia

coeffsB1 = scipy.vstack((B0, B11, B2)) / gain
coeffsB2 = scipy.vstack((B0, B12, B2))
coeffsB3 = scipy.vstack((B0, B13, B2))
coeffsB4 = scipy.vstack((B0, B14, B2))

coeffsA = scipy.vstack((A0, A1, A2))

npoints = 10000
w = scipy.arange(-scipy.pi, scipy.pi, 2*scipy.pi/npoints)

d1 = loudia.freqz(coeffsB1, coeffsA, w)
d2 = loudia.freqz(coeffsB2, coeffsA, w)
d3 = loudia.freqz(coeffsB3, coeffsA, w)
d4 = loudia.freqz(coeffsB4, coeffsA, w)

d = d1 * d2 * d3 * d4

# Plot the frequency response
import pylab

subplots = 2 if plotAngle else 1;


pylab.subplot(subplots,1,1)
if plotColor:
    pylab.plot(w[npoints/2:], loudia.magToDb(abs(d[npoints/2:,:])))
else:
    pylab.plot(w[npoints/2:], loudia.magToDb(abs(d[npoints/2:,:])), c = 'black')

pylab.hold(True)

suma = abs(d[npoints/2:,:]).sum(axis=1)
suma.resize((suma.shape[0], 1))
pylab.plot(w[npoints/2:], loudia.magToDb(suma), c = 'red')

ax = pylab.gca()

# Show half of the spectrum
ax.set_xlim([0, scipy.pi])

# Set the ticks units to radians per second
ticks = ax.get_xticks()
ax.set_xticklabels(['%.2f' % (float(tick) / (scipy.pi * 2.0)) for tick in ticks])

# Set the title and labels
pylab.title('Magnitude of the Frequency Response of a \n Gammatone Filterbank implementation')
pylab.xlabel('Normalized Frequency')
pylab.ylabel('|H(w)| (no unit)')

if plotAngle:
    pylab.subplot(2,1,2)
    if plotColor:
        pylab.plot(w[npoints/2:], scipy.angle(d[npoints/2:,:]))
    else:
        pylab.plot(w[npoints/2:], scipy.angle(d[npoints/2:,:]), c = 'black')
    pylab.title('Angle of the Frequency Response')

    pylab.gca().set_xlim([0, scipy.pi])

pylab.show()

