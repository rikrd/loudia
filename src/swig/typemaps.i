/* -*- C -*-  (not really, but good for syntax highlighting) */

/*
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or   
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

%apply float { Real };

%typecheck(SWIG_TYPECHECK_INTEGER)
	 int, short, long,
 	 unsigned int, unsigned short, unsigned long,
	 signed char, unsigned char,
	 long long, unsigned long long,
	 const int &, const short &, const long &,
 	 const unsigned int &, const unsigned short &, const unsigned long &,
	 const long long &, const unsigned long long &,
	 enum SWIGTYPE,
         bool, const bool & 
{
  $1 = (PyInt_Check($input) || PyLong_Check($input)) ? 1 : 0;
}


%typemap(typecheck,
         prefedence = SWIG_TYPECHECK_FLOAT) 
         Real, 
         const Real,
         Real & {

  $1 = (PyFloat_Check($input) || PyInt_Check($input) || PyLong_Check($input)) ? 1 : 0;

}

%typemap(typecheck,
         prefedence = SWIG_TYPECHECK_FLOAT_ARRAY) 
         MatrixXR, 
         const MatrixXR,
         MatrixXR & {
  $1 = type_match(array_type($input), PyArray_FLOAT);
}

%typemap(typecheck,
         prefedence = SWIG_TYPECHECK_FLOAT_ARRAY) 
         MatrixXC,
         const MatrixXC,
         MatrixXC & {
  $1 = type_match(array_type($input), PyArray_CFLOAT);
}

%typemap(in,
         fragment="NumPy_Fragments") 
         MatrixXR {

    int newObject;
    PyArrayObject * in_array = obj_to_array_contiguous_allow_conversion($input, PyArray_FLOAT, &newObject);

    if( in_array == NULL ){
      PyErr_SetString(PyExc_ValueError,
                      "array must be of type float (dtype = 'float32')");
      
      return NULL;
    }
    
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
    Eigen::Map<MatrixXRscipy> in_matrix(in_data, in_rows, in_cols);

    $1.set(in_matrix);
}

%typemap(in,
         fragment="NumPy_Fragments")
         MatrixXC {

    int newObject;
    PyArrayObject * in_array = obj_to_array_contiguous_allow_conversion($input, PyArray_CFLOAT, &newObject);
    
    if( in_array == NULL ){
      PyErr_SetString(PyExc_ValueError,
                      "array must be of type complex float (dtype = 'complex64')");
      
      return NULL;
    }

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

    Eigen::Map<MatrixXCscipy> in_matrix(in_data, in_rows, in_cols);

    $1.set(in_matrix);
}

%typemap(in, numinputs = 0) 
         MatrixXR* (MatrixXR temp) {

  $1 = &temp;

}

%typemap(argout) 
         MatrixXR* {

  // prepare resulting array
  int dims[] = {(*$1).rows(), (*$1).cols()};
  PyObject * out_array = PyArray_FromDims(2, dims, PyArray_FLOAT);

  if (out_array == NULL){
    PyErr_SetString(PyExc_ValueError,
                    "Unable to create the output array.");
    
    return NULL;
  }
  
  Real* out_data = (Real*)array_data(out_array);
  Eigen::Map<MatrixXRscipy>(out_data, dims[0], dims[1]) = (*$1);

  $result = SWIG_Python_AppendOutput($result, out_array);
}

%typemap(in, numinputs = 0) 
         MatrixXC* (MatrixXC temp) {

  $1 = &temp;

}

%typemap(argout) 
         MatrixXC* {

  // prepare resulting array
  int dims[] = {(*$1).rows(), (*$1).cols()};
  PyObject * out_array = PyArray_FromDims(2, dims, PyArray_CFLOAT);

  if (out_array == NULL){
    PyErr_SetString(PyExc_ValueError,
                    "Unable to create the output array.");
    
    return NULL;
  }
  
  Complex* out_data = (Complex*)array_data(out_array);
  Eigen::Map<MatrixXCscipy>(out_data, dims[0], dims[1]) = (*$1);

  $result = SWIG_Python_AppendOutput($result, out_array);
}
