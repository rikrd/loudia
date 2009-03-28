#!/usr/bin/env python

# Create input
import scipy.signal
import ricaudio

frameSize = 256

a_ones = scipy.ones( frameSize )

# Ricaudio's solution # --------------------------------- #
m_hanning = ricaudio.Window( frameSize, 1 )
m_hamming = ricaudio.Window( frameSize, 3 )

r_hanning = m_hanning.process( a_ones )
r_hamming = m_hamming.process( a_ones )
# -------------------------------------------------------- #

# Scipy solution # ---------------------------------- #
s_hanning = scipy.signal.hanning( frameSize )
s_hamming = scipy.signal.hamming( frameSize )
# -------------------------------------------------------- #

print 'Hann match: ', scipy.allclose( r_hanning, s_hanning )

print 'Hamming match: ',scipy.allclose( r_hamming, s_hamming )

print ( r_hamming - s_hamming ).sum()

#import pylab
#pylab.hold(True)
#pylab.plot(r_hamming.T)
#pylab.plot(s_hamming)
#pylab.show()
