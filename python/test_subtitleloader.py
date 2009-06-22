#!/usr/bin/env python

# Create input
import scipy
import loudia

sampleRate = 120
frameSize = 100
hopSize = 3
loadDuration = 10

# Loudia's solution # --------------------------------- #
loader = loudia.SubtitleLoader("test.srt", sampleRate, frameSize, hopSize, loadDuration)

frames = []

while not loader.isFinished():
    frame = loader.process()

    frames += frame

print frames
