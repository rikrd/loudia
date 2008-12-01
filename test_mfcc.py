#!/usr/bin/env python

# Create input
import scipy
lowFreq = 300.0
highFreq = 16000.0
nBands = 40
samplerate = 44100
spectralLength = 1024
nCoeffs = 13

a1 = scipy.array(scipy.zeros((512, 1)), dtype='f4')
a2 = scipy.array(scipy.ones((512, 1)), dtype='f4')

# CRicaudio's solution # --------------------------------- #
import pycricaudio
m = pycricaudio.MFCC(lowFreq, highFreq, nBands, samplerate, spectralLength, nCoeffs)

b1 = m.process(a1)
b2 = m.process(a2)
# -------------------------------------------------------- #

# Ricaudio's solution # ---------------------------------- #
from aubio import aubioclass
def mfcc_compute(samples):
    mvec = aubioclass.fvec(spectralLength,
                           1)
    
    for ind in range(samples.shape[0]):
        aubioclass.fvec_write_sample(mvec(), float(samples[ind]), 0, ind)
    
    mopick = aubioclass.onsetpick(spectralLength,
                                  spectralLength,
                                  1,
                                  mvec,
                                  0.6,
                                  mode = 'dual',
                                  derivate = False,
                                  dcthreshold = 0)
    
    mfcc = aubioclass.new_aubio_mfcc(spectralLength,
                                     samplerate,
                                     nBands,
                                     nCoeffs,
                                     lowFreq,
                                     highFreq,
                                     1)
    
    mfcc_coeffs = aubioclass.new_fvec(nCoeffs,
                                      1)
    
    (onset, value) = mopick.do(mvec)
    aubioclass.aubio_mfcc_do(mfcc, mopick.myfft(), mfcc_coeffs)
        
    mfcc_coefficients = scipy.zeros((nCoeffs, 1))
    for ind in range(nCoeffs):
        mfcc_coefficients[ind, 0] = aubioclass.fvec_read_sample(mfcc_coeffs, 0, ind)

    aubioclass.del_aubio_mfcc(mfcc)
    del mfcc
    
    aubioclass.del_fvec(mfcc_coeffs)
    del mfcc_coeffs

    return mfcc_coefficients

c1 = mfcc_compute(a1)
c2 = mfcc_compute(a2)
# -------------------------------------------------------- #

print b1
print b2

print c1
print c2

print scipy.allclose(b1,c1)
print scipy.allclose(b2,c2)
