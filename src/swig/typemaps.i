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

%apply float { Real };

%exception {
  try {
    $action
  } catch (LoudiaException e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return NULL;
  }
}

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


%typecheck(SWIG_TYPECHECK_FLOAT) 
           Real,
           const Real,
           Real & {

  $1 = (PyFloat_Check($input) || PyInt_Check($input) || PyLong_Check($input)) ? 1 : 0;

}

%typecheck(SWIG_TYPECHECK_FLOAT_ARRAY) 
         MatrixXR, 
         MatrixXR *,
         const MatrixXR,
         MatrixXR &,
         const MatrixXR & {
  $1 = (array_type($input) == PyArray_FLOAT) || (array_type($input) == PyArray_DOUBLE);
}

%typecheck(SWIG_TYPECHECK_FLOAT_ARRAY) 
         MatrixXC,
         MatrixXC *,
         const MatrixXC,
         MatrixXC &,
         const MatrixXC & {
  $1 = (array_type($input) == PyArray_CFLOAT) || (array_type($input) == PyArray_CDOUBLE);
}

%typemap(in,
         fragment="NumPy_Fragments") 
         const MatrixXR & (MatrixXR temp) {

    // create array from input
    int newObject;
    PyArrayObject * in_array;
    
    switch ( array_type($input) ) {

    case PyArray_LONG:
    case PyArray_DOUBLE:
      in_array = obj_to_array_contiguous_allow_conversion($input, PyArray_DOUBLE, &newObject);
      break;

    case PyArray_INT:
    case PyArray_FLOAT:
      in_array = obj_to_array_contiguous_allow_conversion($input, PyArray_FLOAT, &newObject);
      break;

    default:
      PyErr_SetString(PyExc_ValueError,
                      "array must be of type int, float, long or double");
      
      return NULL;
    }

    if( in_array == NULL ){
      PyErr_SetString(PyExc_ValueError,
                      "array could not be created");
      
      return NULL;
    }
    
    // require one or two dimensions
    npy_intp dims[] = {1, 2};
    require_dimensions_n(in_array, dims, 2);

    // get the dimensions
    int in_rows;
    int in_cols;
    if(array_numdims(in_array) == 2){

      in_rows = array_size(in_array, 0);
      in_cols = array_size(in_array, 1);

    }else{

      in_rows = 1;
      in_cols = array_size(in_array, 0);

    }

    $1 = &temp;

    // prepare the input array
    switch( array_type($input) ) {

    case PyArray_LONG:
    case PyArray_DOUBLE:
      (*$1) = Eigen::Map<MatrixXdscipy>((double*)array_data( in_array ), in_rows, in_cols).cast<Real>();
      break;

    case PyArray_INT:
    case PyArray_FLOAT:
      (*$1) = Eigen::Map<MatrixXfscipy>((float*)array_data( in_array ), in_rows, in_cols).cast<Real>();
      break;
      
    default:
      PyErr_SetString(PyExc_ValueError,
                      "array must be of type int, float, long or double");
      return NULL;
    }
}


%typemap(in,
         fragment="NumPy_Fragments") 
         const MatrixXC & (MatrixXC temp) {

    // create array from input
    int newObject;
    PyArrayObject * in_array;

    switch ( array_type($input) ) {

    case PyArray_LONG:
    case PyArray_DOUBLE:
    case PyArray_CDOUBLE:
      in_array = obj_to_array_contiguous_allow_conversion($input, PyArray_CDOUBLE, &newObject);
      break;

    case PyArray_INT:
    case PyArray_FLOAT:
    case PyArray_CFLOAT:
      in_array = obj_to_array_contiguous_allow_conversion($input, PyArray_CFLOAT, &newObject);
      break;
    
    default:
      PyErr_SetString(PyExc_ValueError,
                      "array must be of type complex int, complex float, complex long or complex double");
      return NULL;
      
    }

    if( in_array == NULL ){
      PyErr_SetString(PyExc_ValueError,
                      "array could not be created");
      
      return NULL;
    }
    
    // require one or two dimensions
    npy_intp dims[] = {1, 2};
    require_dimensions_n(in_array, dims, 2);

    // get the dimensions
    int in_rows;
    int in_cols;
    if(array_numdims(in_array) == 2){

      in_rows = array_size(in_array, 0);
      in_cols = array_size(in_array, 1);

    }else{

      in_rows = 1;
      in_cols = array_size(in_array, 0);

    }

    $1 = &temp;

    // prepare the input array
    switch( array_type($input) ) {

    case PyArray_LONG:
    case PyArray_DOUBLE:
    case PyArray_CDOUBLE:
      (*$1) = Eigen::Map<MatrixXcdscipy>((std::complex<double> *) array_data( in_array ), in_rows, in_cols).cast<Complex>();
      break;

    case PyArray_INT:
    case PyArray_FLOAT:
    case PyArray_CFLOAT:
      (*$1) = Eigen::Map<MatrixXcfscipy>((std::complex<float> *) array_data( in_array ), in_rows, in_cols).cast<Complex>();
      break;
      
    default:
      PyErr_SetString(PyExc_ValueError,
                      "array must be of type complex int, complex float, complex long or complex double");
      return NULL;
    }
}



%typemap(in, numinputs = 0) 
         MatrixXR* (MatrixXR temp) {

  $1 = &temp;

}

%typemap(argout) 
         MatrixXR* {

  // prepare resulting array
  npy_intp dims[] = {$1->rows(), $1->cols()};
  PyObject * out_array = PyArray_SimpleNew(2, dims, PyArray_FLOAT);

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
         MatrixXI* (MatrixXI temp) {

  $1 = &temp;

}

%typemap(argout) 
         MatrixXI* {

  // prepare resulting array
  npy_intp dims[] = {$1->rows(), $1->cols()};
  PyObject * out_array = PyArray_SimpleNew(2, dims, PyArray_INT);

  if (out_array == NULL){
    PyErr_SetString(PyExc_ValueError,
                    "Unable to create the output array.");
    
    return NULL;
  }
  
  Integer* out_data = (Integer*)array_data(out_array);
  Eigen::Map<MatrixXIscipy>(out_data, dims[0], dims[1]) = (*$1);

  $result = SWIG_Python_AppendOutput($result, out_array);
}

%typemap(in, numinputs = 0) 
         MatrixXC* (MatrixXC temp) {

  $1 = &temp;

}

%typemap(argout) 
         MatrixXC* {

  // prepare resulting array
  npy_intp dims[] = {$1->rows(), $1->cols()};
  PyObject * out_array = PyArray_SimpleNew(2, dims, PyArray_CFLOAT);

  if (out_array == NULL){
    PyErr_SetString(PyExc_ValueError,
                    "Unable to create the output array.");
    
    return NULL;
  }
  
  Complex* out_data = (Complex*)array_data(out_array);
  Eigen::Map<MatrixXCscipy>(out_data, dims[0], dims[1]) = (*$1);

  $result = SWIG_Python_AppendOutput($result, out_array);
}

%typemap(in, numinputs = 0) 
         Real* (Real temp) {

  $1 = &temp;

}

%typemap(argout) 
         Real* {

  $result = SWIG_Python_AppendOutput($result, Py_BuildValue("f", *$1));
}


%typemap(in, numinputs = 0) 
         int* (int temp) {

  $1 = &temp;

}

%typemap(argout) 
         int* {

  $result = SWIG_Python_AppendOutput($result, Py_BuildValue("i", *$1));
}
