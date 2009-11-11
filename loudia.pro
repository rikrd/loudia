CONFIG -= qt \
    gui \
    core
TEMPLATE = lib
QMAKE_CXX = ccache g++
QMAKE_CC = ccache gcc
HEADERS += src/Window.h \
    src/VoiceActivityDetection.h \
    src/Utils.h \
    src/Unwrap.h \
    src/Typedefs.h \
    src/SpectralWhitening.h \
    src/SpectralReassignment.h \
    src/SpectralODFPhase.h \
    src/SpectralODFMKL.h \
    src/SpectralODFHFC.h \
    src/SpectralODFFlux.h \
    src/SpectralODFComplex.h \
    src/SpectralODFCOG.h \
    src/SpectralODFBase.h \
    src/SpectralODF.h \
    src/SpectralNoiseSuppression.h \
    src/Resample.h \
    src/PitchSaliency.h \
    src/PitchInverseProblem.h \
    src/PitchACF.h \
    src/PeakTracking.h \
    src/PeakSynthesize.h \
    src/PeakInterpolationComplex.h \
    src/PeakInterpolation.h \
    src/PeakDetectionComplex.h \
    src/PeakDetection.h \
    src/PeakCOG.h \
    src/OnsetComplex.h \
    src/NMF.h \
    src/MFCC.h \
    src/MelScales.h \
    src/MelBands.h \
    src/Meddis.h \
    src/MatrixBaseAddons.h \
    src/LPCResidual.h \
    src/LPC.h \
    src/INMF.h \
    src/IFFTComplex.h \
    src/IFFT.h \
    src/FunctorsAddons.h \
    src/FrameCutter.h \
    src/FilterUtils.h \
    src/Filter.h \
    src/FFTComplex.h \
    src/FFT.h \
    src/Debug.h \
    src/DCT.h \
    src/CwiseAddons.h \
    src/Correlation.h \
    src/BarkBands.h \
    src/Bands.h \
    src/BandFilter.h \
    src/Autocorrelation.h \
    src/AudioLoader.h \
    src/AOK.h
SOURCES += src/Window.cpp \
    src/VoiceActivityDetection.cpp \
    src/Utils.cpp \
    src/Unwrap.cpp \
    src/SpectralWhitening.cpp \
    src/SpectralReassignment.cpp \
    src/SpectralODFPhase.cpp \
    src/SpectralODFMKL.cpp \
    src/SpectralODFHFC.cpp \
    src/SpectralODFFlux.cpp \
    src/SpectralODFComplex.cpp \
    src/SpectralODFCOG.cpp \
    src/SpectralODFBase.cpp \
    src/SpectralODF.cpp \
    src/SpectralNoiseSuppression.cpp \
    src/Resample.cpp \
    src/PitchSaliency.cpp \
    src/PitchInverseProblem.cpp \
    src/PitchACF.cpp \
    src/PeakTracking.cpp \
    src/PeakSynthesize.cpp \
    src/PeakInterpolationComplex.cpp \
    src/PeakInterpolation.cpp \
    src/PeakDetectionComplex.cpp \
    src/PeakDetection.cpp \
    src/PeakCOG.cpp \
    src/OnsetComplex.cpp \
    src/NMF.cpp \
    src/MFCC.cpp \
    src/MelScales.cpp \
    src/MelBands.cpp \
    src/Meddis.cpp \
    src/LPCResidual.cpp \
    src/LPC.cpp \
    src/INMF.cpp \
    src/IFFTComplex.cpp \
    src/IFFT.cpp \
    src/FrameCutter.cpp \
    src/FilterUtils.cpp \
    src/Filter.cpp \
    src/FFTComplex.cpp \
    src/FFT.cpp \
    src/DCT.cpp \
    src/Correlation.cpp \
    src/BarkBands.cpp \
    src/Bands.cpp \
    src/BandFilter.cpp \
    src/Autocorrelation.cpp \
    src/AudioLoader.cpp \
    src/AOK.cpp
LIBS += -lfftw3f \
        -lavcodec \
        -lavutil \
        -lavformat \
        -lsamplerate
INCLUDEPATH += ../eigen2 \
    src
