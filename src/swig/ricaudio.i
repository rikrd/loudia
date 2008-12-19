%module ricaudio

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
import_array();
%}

%{
#include "melbands.h"
%}


#include "melbands.h"

