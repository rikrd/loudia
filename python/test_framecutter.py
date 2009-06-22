#!/usr/bin/env python

# Create input
import scipy
import loudia

inputSize = 12
frameSize = 5
hopSize = 3

# Loudia's solution # --------------------------------- #
window = loudia.FrameCutter(inputSize, frameSize, hopSize, 5)

inputFrames = 5
streamFrames = scipy.arange(inputFrames*inputSize).reshape((inputFrames, inputSize)).T

print streamFrames

frames = []
for i in range(inputFrames):
    streamFrame = streamFrames[:, i]
    streamFrame = streamFrame.reshape((inputSize, 1))
    produced, newFrames = window.process(streamFrame)

    frames += list(newFrames[:produced, :])

print frames
