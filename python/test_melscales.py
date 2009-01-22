#!/usr/bin/env python

import scipy
import ricaudio

lowFreq = 133.0
highFreq = 22050.0

f_m = lowFreq
l = ricaudio.linearToMelGreenwood1990( f_m )
f_c = ricaudio.melToLinearGreenwood1990( l )

print 'Low Freq Good:', f_m, ' Hz'
print 'Low Mel Good:', l, ' mels' 

print scipy.allclose(f_c, f_m)

f_m = [lowFreq]
l = ricaudio.linearToMelMatrixGreenwood1990( f_m )
f_c = ricaudio.melToLinearMatrixGreenwood1990( l )


print 'Low Freq Bad:', f_m, 'Hz'
print 'Low Mel Bad:', l, 'mels' 
print 'Low Freq Calculated Bad', f_c, 'Hz'

print scipy.allclose(f_c, f_m)


f_m = highFreq
l = ricaudio.linearToMelGreenwood1990( f_m )
f_c = ricaudio.melToLinearGreenwood1990( l )

print 'High Freq Good:', f_m, 'Hz'
print 'High Mel Good:', l, 'mels' 

print scipy.allclose(f_c, f_m)

f_m = [highFreq]
l = ricaudio.linearToMelMatrixGreenwood1990( f_m )
f_c = ricaudio.melToLinearMatrixGreenwood1990( l )


print 'High Freq Bad:', f_m, 'Hz'
print 'High Mel Bad:', l, 'mels'
print 'High Freq Calculated Bad', f_c, 'Hz'

print scipy.allclose(f_c, f_m)

