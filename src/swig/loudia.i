/* -*- C -*-  (not really, but good for syntax highlighting) */

/*
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 3 of the License, or   
** (at your option) any later version.                                 
**                                                                     
** This program is distributed in the hope that it will be useful,     
** but WITHOUT ANY WARRANTY; without even the implied warranty of      
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
** GNU General Public License for more details.                        
**                                                                     
** You should have received a copy of the GNU General Public License   
** along with this program; if not, write to the Free Software         
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
*/

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
#include "BarkBands.h"
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
#include "VoiceActivityDetection.h"

#include "MelScales.h"
#include "Utils.h"
#include "FilterUtils.h"
%}

%include "typemaps.i"

%include "Filter.h"
%include "DCT.h"
%include "Window.h"
%include "MelBands.h"
%include "BarkBands.h"
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
%include "VoiceActivityDetection.h"

%include "MelScales.h"
%include "Utils.h"
%include "FilterUtils.h"

%module loudia
