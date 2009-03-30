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
#include "PeakDetection.h"
#include "PeakDetectionComplex.h"
#include "PeakCOG.h"
#include "PeakInterpolation.h"
#include "PeakInterpolationComplex.h"
#include "PeakTracking.h"
#include "PeakSynthesize.h"
#include "BandFilter.h"
#include "Unwrap.h"
#include "LPC.h"
#include "LPCResidual.h"
#include "NMF.h"
#include "INMF.h"
#include "Resample.h"
#include "Correlation.h"
#include "Autocorrelation.h"
#include "SpectralNoiseSuppression.h"
#include "SpectralWhitening.h"
#include "SpectralODF.h"
#include "PitchSaliency.h"
#include "PitchACF.h"
#include "PitchInverseProblem.h"

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
%include "PeakDetection.h"
%include "PeakDetectionComplex.h"
%include "PeakCOG.h"
%include "PeakInterpolation.h"
%include "PeakInterpolationComplex.h"
%include "PeakTracking.h"
%include "PeakSynthesize.h"
%include "BandFilter.h"
%include "Unwrap.h"
%include "LPC.h"
%include "LPCResidual.h"
%include "NMF.h"
%include "INMF.h"
%include "Resample.h"
%include "Correlation.h"
%include "Autocorrelation.h"
%include "SpectralNoiseSuppression.h"
%include "SpectralWhitening.h"
%include "SpectralODF.h"
%include "PitchSaliency.h"
%include "PitchACF.h"
%include "PitchInverseProblem.h"

%include "MelScales.h"
%include "Utils.h"
%include "FilterUtils.h"

%module loudia
