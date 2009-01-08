#!/usr/bin/env python

import ricaudio
import scipy

a = [[1,  2,  6,  4],
     [1,  2,  6,  4],
     [10,  20,  60,  40]]

d = ricaudio.ODFComplex(4, 16, ricaudio.Window.RECTANGULAR, False)

print 'start ricaudio'
#for i in range(10000):
b1 = d.process(a)
print 'end'


def princarg(phase):
    scipy.add(phase, scipy.pi, phase)
    scipy.remainder(phase,-2.0*scipy.pi, phase)
    scipy.add(phase, scipy.pi, phase)
    return phase

predictMethod = 'const_lineal'
if predictMethod == 'const_lineal':
    def predict(passedMags, passedPhases):
        predictMagnitude = passedMags[-2]
        predictPhase = passedPhases[-2] + ( passedPhases[-2] - passedPhases[-3] )
        return (predictMagnitude, predictPhase)

elif predictMethod == 'lineal_lineal':
    def predict(passedMags, passedPhases):
        predictMagnitude = (passedMags[-2] + ( passedMags[-2] - passedMags[-3] ))
        predictPhase = passedPhases[-2] + ( passedPhases[-2] - passedPhases[-3] )
        return (predictMagnitude, predictPhase)
else:
    raise RuntimeError, "The predict method selected is not implemented.  Please select one of the following: 'const_lineal' or 'lineal_lineal'."

distanceMethod = 'hypot'
if distanceMethod == 'euclidean':
    def distance(predictFft, currentFft):
        return scipy.hypot(currentFft.real - predictFft.real,
                           currentFft.imag - predictFft.imag)
    
elif distanceMethod == 'euclidean_weighted':
    def distance(predictFft, currentFft):
        return scipy.hypot(currentFft.real - predictFft.real,
                           currentFft.imag - predictFft.imag) * scipy.absolute(currentFft)
    
elif distanceMethod == 'hypot':
    def distance(predictFft, currentFft):
        return (scipy.hypot(predictFft.real, currentFft.real) +
                scipy.hypot(predictFft.imag, currentFft.imag))
    
else:
    raise RuntimeError, "The distance method selected is not implemented.  Please select one of the following: 'euclidean', 'euclidean_weighted' or 'hypot'."


j = scipy.complex64(1.0j)

passedPhases = scipy.zeros((1,1))
passedMags = scipy.zeros((1,1))
diff = scipy.zeros((1,1))
upsteps = scipy.zeros((1,1))
downsteps = scipy.zeros((1,1))
shift = scipy.zeros((1,1))

weights = scipy.ones((1,1))


# Fast angle that doesn't check for errors
def fastangle(z):
    # returns the angle of the complex numbers in array z
    return scipy.arctan2(z.imag, z.real)

# Use a different unwarp since scipy's is 2x slower
def fastunwrap(thetaArray, discont = scipy.pi):
    # takes an array of theta values
    # returns the data in unwrapped form (unwrapping over the axis == 1)
    diff = scipy.zeros_like(thetaArray)
    diff[1:,:] = scipy.diff(thetaArray, axis = 0)
    upSteps = diff > discont
    downSteps = diff < -discont
    shift = scipy.cumsum(upSteps, axis = 0) - scipy.cumsum(downSteps, axis = 0)
    return thetaArray - 2.0*discont*shift

# Fast diff uses ufuncs
def veryfastdiff(a, output):
    # returns the 1st order difference along the axis
    scipy.add(-a[:-1], a[1:], output)

# Use a different unwarp since scipy's is 2x slower
def veryfastunwrap(thetaArray, diff, upsteps, downsteps, shift, output, discont = scipy.pi):
    # takes an array of theta values
    # returns the data in unwrapped form (unwrapping over the axis == 1)
    diff[::] = 0
    veryfastdiff(thetaArray, diff[1:,:])

    scipy.greater(diff, discont, upsteps)
    scipy.cumsum(upsteps, axis = 0, out = upsteps)

    scipy.less(diff, -discont, downsteps)
    scipy.cumsum(downsteps, axis = 0, out = downsteps)

    scipy.subtract(upsteps, downsteps, shift)

    scipy.add(thetaArray, - 2.0*discont*shift, output)

def wmean(vector, weights):
    return (vector*weights).sum()/weights.sum()

def veryfastSpectralPredictionError(vector, cursor):
    passedFft = vector[:,:vector.shape[1]/2]
    
    # Initialize buffers
    if passedPhases.shape != passedFft.shape:
        passedMags.resize(passedFft.shape)
        passedPhases.resize(passedFft.shape)
        diff.resize(passedFft.shape)
        upsteps.resize(passedFft.shape)
        downsteps.resize(passedFft.shape)
        shift.resize(passedFft.shape)
    
    passedWrappedPhases = fastangle(passedFft)
    
    veryfastunwrap(passedWrappedPhases, diff, upsteps, downsteps, shift, passedPhases)
    
    scipy.absolute(passedFft, passedMags)
    
    # Calculate the predicted complex fft
    predictMagnitude, predictPhase = predict(passedMags, passedPhases)
    predictFft = predictMagnitude * scipy.exp( predictPhase * j )
    
    # Calculate prediction error
    predictionError = distance(predictFft, passedFft[-1])
       
    # Integrate the prediction error of all bins
    predictionError = predictionError[1:].mean()
    
    return predictionError

def spectralPredictionError(vector, cursor):

    passedFft = vector[:,:vector.shape[1]/2]
    
    # Unwrap the phases
    passedPhases = fastunwrap(fastangle(passedFft))
    passedMags = scipy.absolute(passedFft)
    
    #print passedFft
    #print passedMags
    #print passedPhases
    
    # Calculate the predicted complex fft
    predictMagnitude, predictPhase = predict(passedMags, passedPhases)
    predictFft = predictMagnitude * scipy.exp( predictPhase * j )
   
    # Calculate prediction error
    predictionError = distance(predictFft, passedFft[-1])
        
    # The energy should have the weight of the other bins
    weights = scipy.ones(predictionError.shape)
    weights[0] = 0.0
    
    # Integrate the prediction error of all bins
    predictionError = wmean(predictionError, weights)
        
    return predictionError

print 'start scipy'
#for i in range(10000):
ffted = scipy.fft(a, 16)
b2 = spectralPredictionError(ffted, 0)
print 'end'

print b1[0,0]
print b2
