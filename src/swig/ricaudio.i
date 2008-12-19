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
%}

%typemap(in,
         fragment="NumPy_Fragments") MatrixXR {
    int newObject;
    PyArrayObject * in_array = obj_to_array_contiguous_allow_conversion($input, PyArray_FLOAT, &newObject);
    
    int dims[] = {1, 2};
    require_dimensions_n(in_array, dims, 2);

    int in_rows;
    int in_cols;

    if(array_numdims(in_array) == 2){

      in_rows = array_size(in_array, 0);
      in_cols = array_size(in_array, 1);

    }else{

      in_rows = 1;
      in_cols = array_size(in_array, 0);

    }

    // prepare the input array  
    Real* in_data = (Real*)array_data(in_array);

    MatrixXR in_matrix = Eigen::Map<MatrixXRscipy>(in_data, in_rows, in_cols);

    $1 = in_matrix;
}

%include "filter.h"

%module ricaudio
