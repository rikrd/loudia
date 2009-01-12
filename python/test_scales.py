#!/usr/bin/env python

import scipy
import ricaudio

lowFreq = 33.0
highFreq = 22000.0
nBands = 34
samplerate = 44100.0
spectralLength = 2048

m = ricaudio.MelBands(lowFreq, highFreq, nBands, samplerate, spectralLength)

f_m = 33.3
l = m.linearToMelGreenwood1990( f_m )
f_c = m.melToLinearGreenwood1990( l )

print scipy.allclose(f_c, f_m)

f_m = [33.3]
l = m.linearToMelMatrixGreenwood1990( f_m )
f_c = m.melToLinearMatrixGreenwood1990( l )

print scipy.allclose(f_c, f_m)

