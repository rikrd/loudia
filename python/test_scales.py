#!/usr/bin/env python

import scipy
import ricaudio

lowFreq = 133.0
highFreq = 22050.0
nBands = 34
samplerate = 44100.0
spectralLength = 2048

m = ricaudio.MelBands(lowFreq, highFreq, nBands, samplerate, spectralLength)

f_m = lowFreq
l = m.linearToMelGreenwood1990( f_m )
f_c = m.melToLinearGreenwood1990( l )

print 'Low Freq Good:', f_m, ' Hz'
print 'Low Mel Good:', l, ' mels' 

print scipy.allclose(f_c, f_m)

f_m = [lowFreq]
l = m.linearToMelMatrixGreenwood1990( f_m )
f_c = m.melToLinearMatrixGreenwood1990( l )


print 'Low Freq Bad:', f_m, 'Hz'
print 'Low Mel Bad:', l, 'mels' 
print 'Low Freq Calculated Bad', f_c, 'Hz'

print scipy.allclose(f_c, f_m)


f_m = highFreq
l = m.linearToMelGreenwood1990( f_m )
f_c = m.melToLinearGreenwood1990( l )

print 'High Freq Good:', f_m, 'Hz'
print 'High Mel Good:', l, 'mels' 

print scipy.allclose(f_c, f_m)

f_m = [highFreq]
l = m.linearToMelMatrixGreenwood1990( f_m )
f_c = m.melToLinearMatrixGreenwood1990( l )


print 'High Freq Bad:', f_m, 'Hz'
print 'High Mel Bad:', l, 'mels'
print 'High Freq Calculated Bad', f_c, 'Hz'

print scipy.allclose(f_c, f_m)

