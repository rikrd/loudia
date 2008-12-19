%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
import_array();
%}

%module ricaudio

%typemap(in) MatrixXR {
    int newObject;
    PyArrayObject * in_array = obj_to_array_contiguous_allow_conversion($input, PyArray_FLOAT, &newObject);
    
    // check that the input array has 2 dimensions
    if (in_array->nd != 2 || in_array->descr->type_num != PyArray_FLOAT) { 
      PyErr_SetString(PyExc_ValueError,
                      "array must be two-dimensional and of type float");
    
      return NULL;
    }

    // prepare the input array
    int in_samples = in_array->dimensions[0];
    int in_channels = in_array->dimensions[1];
  
    int in_rows = in_array->dimensions[0];
    int in_cols = in_array->dimensions[1];
    Real* in_data = (Real*)array_data(in_array);

    MatrixXR in_matrix = Eigen::Map<MatrixXRscipy>(in_data, in_rows, in_cols);

    $1 = in_matrix;
}


%{
#include "filter.h"
%}


%include "filter.h"

