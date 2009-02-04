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

#include <stdio.h>

#include "filter.h"
#include "dct.h"
#include "window.h"
#include "melbands.h"
#include "bands.h"
#include "fft.h"
#include "fftcomplex.h"
#include "ifft.h"
#include "mfcc.h"
#include "aok.h"
#include "meddis.h"
#include "spectralreassignment.h"
#include "peakpick.h"
#include "peakdetect.h"
#include "peakcog.h"
#include "peakinterpolate.h"
#include "peakcontinue.h"
#include "peaksynthesize.h"
#include "unwrap.h"
#include "odfcomplex.h"
#include "odf.h"
#include "odfcog.h"
#include "onsetcomplex.h"
#include "lpc.h"
#include "nmf.h"

#include "melscales.h"
#include "utils.h"

using namespace std;
%}

%include "typemaps.i"

%include "filter.h"
%include "dct.h"
%include "window.h"
%include "bands.h"
%include "melbands.h"
%include "fft.h"
%include "fftcomplex.h"
%include "ifft.h"
%include "mfcc.h"
%include "aok.h"
%include "meddis.h"
%include "spectralreassignment.h"
%include "peakpick.h"
%include "peakdetect.h"
%include "peakcog.h"
%include "peakinterpolate.h"
%include "peakcontinue.h"
%include "peaksynthesize.h"
%include "unwrap.h"
%include "odfcomplex.h"
%include "odf.h"
%include "odfcog.h"
%include "onsetcomplex.h"
%include "lpc.h"
%include "nmf.h"

%include "melscales.h"
%include "utils.h"

%module ricaudio
