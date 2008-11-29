#!/usr/bin/env python

# Create input
import scipy
samplerate = 44100
channels = 30
a1 = scipy.array(scipy.random.random((4096, channels)), dtype='f4')
a2 = scipy.array(scipy.random.random((4096, channels)), dtype='f4')
substractSpont = True

# CRicaudio's solution # --------------------------------- #
import pycricaudio
m = pycricaudio.Meddis(samplerate, channels)

b1 = m.process(a1)
b2 = m.process(a2)
# -------------------------------------------------------- #

# Ricaudio's solution # ---------------------------------- #
M = 1.
A = 5.
B = 300.
g = 2000.
y = 5.05
l = 2500.
r = 6580.
x = 66.31
h = 50000.

# Internal constants
dt = 1/float(samplerate)
gdt = g*dt
ydt = y*dt
ldt = l*dt
rdt = r*dt
xdt = x*dt

# Initial values
init = [True]
kt = [None]
spont = [None]
c = [None]
q = [None]
w = [None]
zeroVector = [None]

#print 'Start meddis...'
def meddis_compute(samples):
    nchannels = samples.shape[1]
    
    if init[0]:
        kt[0] = g*A/(A+B)
        spont[0] = M*y*kt[0]/(l*kt[0]+y*(l+r))
        c[0] = spont[0] * scipy.ones(nchannels)
        q[0] = c[0]*(l+r)/kt[0]
        w[0] = c[0]*r/x
        zeroVector[0] = scipy.zeros(nchannels)
        init[0] = False
        
    def meddis_iteration(row):
        limitedSt = scipy.maximum(row + A, 0.)
        kt[0] = gdt * limitedSt / (limitedSt + B)
        replenish = scipy.maximum(ydt * (M-q[0]), zeroVector[0])
        eject = kt[0] * q[0]
        loss = ldt * c[0]
        reuptake = rdt * c[0]
        reprocess = xdt * w[0]
        
        q[0] += replenish - eject + reprocess
        c[0] += eject - loss - reuptake
        w[0] += reuptake - reprocess
        
        # Now iterate through each time slice of the data.  Use the
        # max function to implement the "if (0>" test.
        out = h * c[0]
        
        if substractSpont:
            out = scipy.maximum(0., out - spont[0])
            
        return out
    
    return scipy.apply_along_axis(meddis_iteration, 1, samples)


c1 = meddis_compute(a1)
c2 = meddis_compute(a2)
# -------------------------------------------------------- #

print scipy.allclose(b1,c1)
print scipy.allclose(b2,c2)
