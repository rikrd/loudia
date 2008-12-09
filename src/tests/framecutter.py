#!/usr/bin/env python

import pylab
import sys

filein = sys.argv[1]
fileout = sys.argv[2]
windowsize = int(sys.argv[3])
framesize = int(2.42*windowsize + 3)
hopsize = 1

ls = ['0']*(framesize-1) + [l.strip().split()[0] for l in open(filein)]

frames = []
end = False
i = 0

while not end:
    newframe = ls[i : i + framesize]

    if len(newframe) == 0:
        break
    
    if len(newframe) < framesize:
        newframe = newframe + ['0.000000']*(framesize - len(newframe))
        
    frames.append(newframe)
    i += hopsize

"""
pylab.ion()
for frame in frames:
    pylab.figure(1)
    pylab.cla()
    pylab.plot([float(a) for a in frame])
    
pylab.ioff()
"""

f = open(fileout, 'w')
for frame in frames:
    f.write(' '.join(frame) + '\n')

f.close()
