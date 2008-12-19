%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
import_array();
%}

%module ricaudio

%{
#include "filter.h"
%}


%include "filter.h"

