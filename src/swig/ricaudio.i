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

#include "filter.h"
#include "dct.h"
#include "window.h"
#include "melbands.h"
#include "bands.h"
#include "fft.h"
#include "mfcc.h"
#include "aok.h"
#include "meddis.h"
#include "spectralreassignment.h"
#include "peakpick.h"
#include "peakdetect.h"
#include "peakinterpolate.h"
#include "chebyshev.h"
#include "unwrap.h"
#include "odfcomplex.h"
#include "onsetcomplex.h"

#include "utils.h"
%}

%include "typemaps.i"

%include "filter.h"
%include "dct.h"
%include "window.h"
%include "bands.h"
%include "melbands.h"
%include "fft.h"
%include "mfcc.h"
%include "aok.h"
%include "meddis.h"
%include "spectralreassignment.h"
%include "peakpick.h"
%include "peakdetect.h"
%include "peakinterpolate.h"
%include "chebyshev.h"
%include "utils.h"
%include "unwrap.h"
%include "odfcomplex.h"
%include "onsetcomplex.h"

%module ricaudio
