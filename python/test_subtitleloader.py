#!/usr/bin/env python

# Create input
import scipy
import loudia

sampleRate = 120.0
frameSize = 100
hopSize = 3
loadDuration = 10.0

# Loudia's solution # --------------------------------- #
loader = loudia.SubtitleLoader()
loader.setFilename("test.srt", False)
loader.setFrameSize(frameSize, False)
loader.setHopSize(hopSize, False)
loader.setLoadDuration(loadDuration, False)
loader.setup()

frames = []
while not loader.isFinished():
    frame = loader.process()
    print frame
    frames += frame

print frames
