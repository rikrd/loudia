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

using namespace std;
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

    $1.set(in_matrix);
}

%typemap(in,
         fragment="NumPy_Fragments") MatrixXC {
    int newObject;
    PyArrayObject * in_array = obj_to_array_contiguous_allow_conversion($input, PyArray_CFLOAT, &newObject);
    
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
    Complex* in_data = (Complex*)array_data(in_array);

    MatrixXC in_matrix = Eigen::Map<MatrixXCscipy>(in_data, in_rows, in_cols);

    $1.set(in_matrix);
}

%typemap(in) MatrixXR* (MatrixXR temp) {
  $1 = &temp;
}

%typemap(argout) MatrixXR* {
    // prepare resulting array
  int dims[] = {(*$1).rows(), (*$1).cols()};
  PyObject * out_array = PyArray_FromDims(2, dims, PyArray_FLOAT);

  if (out_array == NULL){
    PyErr_SetString(PyExc_ValueError,
                    "Unable to create the output array.");
    
    return NULL;
  }

  int out_rows = ((PyArrayObject*)out_array)->dimensions[0];
  int out_cols = ((PyArrayObject*)out_array)->dimensions[1];
  Real* out_data = (Real*)array_data(out_array);

  Eigen::Map<MatrixXRscipy>(out_data, out_rows, out_cols) = (*$1);

  $result = out_array;
}

%include "filter.h"
%include "dct.h"

%module ricaudio
