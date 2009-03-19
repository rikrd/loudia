/* -*- C -*-  (not really, but good for syntax highlighting) */

%{
#define SWIG_FILE_WITH_INIT
%}


%include "numpy.i"

%init %{
import_array();
%}

%{
#include <Eigen/Core>
#include <Eigen/Array>

#include "Filter.h"
#include "DCT.h"
#include "Window.h"
#include "MelBands.h"
#include "Bands.h"
#include "FFT.h"
#include "FFTComplex.h"
#include "IFFT.h"
#include "MFCC.h"
#include "AOK.h"
#include "Meddis.h"
#include "SpectralReassignment.h"
#include "PeakPick.h"
#include "PeakDetect.h"
#include "PeakDetectComplex.h"
#include "PeakCOG.h"
#include "PeakInterpolate.h"
#include "PeakInterpolateComplex.h"
#include "PeakContinue.h"
#include "PeakSynthesize.h"
#include "IIRFilter.h"
#include "LowPass.h"
#include "HighPass.h"
#include "BandPass.h"
#include "BandStop.h"
#include "Unwrap.h"
#include "ODFComplex.h"
#include "ODF.h"
#include "ODFCOG.h"
#include "OnsetComplex.h"
#include "LPC.h"
#include "LPCResidual.h"
#include "NMF.h"
#include "INMF.h"
#include "Resample.h"
#include "Correlation.h"
#include "Autocorrelation.h"
#include "SpectralNoiseSuppression.h"
#include "SpectralWhitening.h"
#include "PitchSaliency.h"
#include "PitchACF.h"

#include "MelScales.h"
#include "Utils.h"
#include "FilterUtils.h"
%}

%include "typemaps.i"

%include "Filter.h"
%include "DCT.h"
%include "Window.h"
%include "MelBands.h"
%include "Bands.h"
%include "FFT.h"
%include "FFTComplex.h"
%include "IFFT.h"
%include "MFCC.h"
%include "AOK.h"
%include "Meddis.h"
%include "SpectralReassignment.h"
%include "PeakPick.h"
%include "PeakDetect.h"
%include "PeakDetectComplex.h"
%include "PeakCOG.h"
%include "PeakInterpolate.h"
%include "PeakInterpolateComplex.h"
%include "PeakContinue.h"
%include "PeakSynthesize.h"
%include "IIRFilter.h"
%include "LowPass.h"
%include "HighPass.h"
%include "BandPass.h"
%include "BandStop.h"
%include "Unwrap.h"
%include "ODFComplex.h"
%include "ODF.h"
%include "ODFCOG.h"
%include "OnsetComplex.h"
%include "LPC.h"
%include "LPCResidual.h"
%include "NMF.h"
%include "INMF.h"
%include "Resample.h"
%include "Correlation.h"
%include "Autocorrelation.h"
%include "SpectralNoiseSuppression.h"
%include "SpectralWhitening.h"
%include "PitchSaliency.h"
%include "PitchACF.h"

%include "MelScales.h"
%include "Utils.h"
%include "FilterUtils.h"

%module ricaudio
