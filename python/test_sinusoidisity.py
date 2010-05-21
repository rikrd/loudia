#!/usr/bin/env python

# Create input
import scipy
import loudia
import pylab
from scikits import audiolab
from viterbi import *

# Synthesize
filename = 'test-sinusoidisity.wav'

sampleRate = 44100.0    # samples / secs

duration = 2.0        # secs

frequency = 440.0      # Hz

vibratoAmplitude = frequency * (pow(2, 1.0/12.0) - 1.0) # Hz
vibratoFrequency = 2    # Hz

sampleCount = duration * sampleRate

noise = (scipy.random.random(sampleCount) - 0.5) * 2

t = scipy.arange(sampleCount) / sampleRate
phi_t = 2.0 * scipy.pi * frequency * t - vibratoAmplitude / vibratoFrequency * scipy.cos(2.0*scipy.pi*vibratoFrequency*t)
sine = scipy.sin(phi_t)

signal = sine + 10*noise

# Writer sinus
audiolab.wavwrite(signal, filename, fs=sampleRate, enc='pcm16')

# Cut into frames
frameSize = 4096                 # samples
hopSize = frameSize/4            # samples
fftSize = frameSize*2            # samples

miniHopSize = frameSize/4        # samples

print miniHopSize

plotSize = 256

frameCount = 0

#windowType = loudia.Window.BLACKMANHARRIS
#windowType = loudia.Window.HANNING
windowType = loudia.Window.HAMMING
#windowType = loudia.Window.RECTANGULAR

window = loudia.Window(frameSize, windowType)
fft = loudia.FFT(fftSize, True)
unwrapper = loudia.Unwrap()

def process(frame):
    spec = fft.process(window.process(frame))
    return scipy.absolute(spec), scipy.angle(spec)


print 'fftSize: ', fftSize


def mapError(error):
    error[error > scipy.pi] = (2*scipy.pi - error[error>scipy.pi])
    return error

errors = []
cursor = 2*miniHopSize
while cursor < (sampleCount - frameSize):
    frameCurrent = signal[cursor:cursor+frameSize]
    framePast = signal[cursor-miniHopSize:cursor-miniHopSize+frameSize]
    framePast2 = signal[cursor-miniHopSize-miniHopSize:cursor-miniHopSize-miniHopSize+frameSize]

    cAbs, cAng = process(frameCurrent)
    pAbs, pAng = process(framePast)
    p2Abs, p2Ang = process(framePast2)

    angs = scipy.vstack((p2Ang, pAng, cAng))

    # Perform the unwrapping of the three frames
    angsUnwrapped = unwrapper.process(angs)

    # Perform the phase prediction and compare to the 3rd frame
    predicted = 2*angsUnwrapped[1,:] - angsUnwrapped[0,:]

    error_linearity = abs(predicted - angsUnwrapped[2,:])

    error_slope = abs((angsUnwrapped[1,:] + 2 * scipy.pi * scipy.arange(angsUnwrapped.shape[1])/(angsUnwrapped.shape[1]-1) * miniHopSize/2.0) - angsUnwrapped[2,:]) % (2*scipy.pi)

    error_linearity = mapError(error_linearity)
    error_slope = mapError(error_slope)

    errors.append(error_linearity + error_slope)

    # The prediction error will be the inverse of the sinusoidity
    if frameCount == 40:
        pylab.figure()
        pylab.subplot(311)
        pylab.plot(loudia.magToDb(cAbs)[0, :plotSize])
        pylab.subplot(312)
        pylab.plot(angsUnwrapped[-1, :plotSize])
        pylab.subplot(313)
        pylab.plot(error_linearity[:plotSize], label='linearity')
        pylab.hold(True)
        pylab.plot(error_slope[:plotSize], label='slope', c='g')

    frameCount += 1

    cursor += hopSize

errors = scipy.array( errors )

import numpy
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]


#errors = errors.max() - errors

smooth_len = 5
errors_smoothed = scipy.zeros_like(errors)
for i in range(errors_smoothed.shape[1]):
    errors_smoothed[:, i] = smooth(errors[:, i], window_len=smooth_len)

for i in range(errors_smoothed.shape[0]):
    errors_smoothed[i, :] = smooth(errors_smoothed[i, :], window_len=smooth_len)


# Perform the ridge detection
def ridgeDetection(image):
    from scikits.image import transform
    out, angles, d = transform.hough(errors)
    pylab.figure()
    pylab.imshow(out, cmap=pylab.cm.bone)
    pylab.xlabel('Angle (degree)')
    pylab.ylabel('Distance %d (pixel)' % d[0])


#ridgeDetection(errors)

print 'Viterbi...'
bestPaths = viterbi_multiple(errors[:, :plotSize])
bestSmoothedPaths = viterbi_multiple(errors_smoothed[:, :plotSize])

print 'Plot...'
pylab.figure()
pylab.subplot(211)
pylab.imshow( errors.T[:plotSize,:], origin='lower', interpolation='nearest')
pylab.hold(True)
pylab.plot( scipy.array(bestPaths).T )
pylab.gca().set_xlim([0, errors.T[:plotSize,:].shape[1]-1])
pylab.gca().set_ylim([0, errors.T[:plotSize,:].shape[0]-1])
pylab.subplot(212)
pylab.imshow( errors_smoothed.T[:plotSize,:], origin='lower', interpolation='nearest')
pylab.hold(True)
pylab.plot( scipy.array(bestSmoothedPaths).T )
pylab.gca().set_xlim([0, errors_smoothed.T[:plotSize,:].shape[1]-1])
pylab.gca().set_ylim([0, errors_smoothed.T[:plotSize,:].shape[0]-1])
pylab.show()

