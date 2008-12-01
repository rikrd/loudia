/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com.com>
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

#include "pyfilter.h"

#include "Eigen/Core"
#include "Eigen/Array"

#define NO_IMPORT_ARRAY
#include "ndarrayobject.h"

#include "typedefs.h"

using namespace std;

PyMethodDef PyFilter_methods[] = {
  { "process",            (PyCFunction)PyFilter::process, METH_VARARGS,
    "Filter.process(input_array) processes the input samples \"input_array\"" },
  { NULL }  /* Sentinel */
};

PyTypeObject PyFilterType = {
  PyObject_HEAD_INIT(NULL)
  0,                         // ob_size
  "cricaudio.Filter",           // tp_name
  sizeof(PyFilter),            // tp_basicsize
  0,                         // tp_itemsize
  PyFilter::dealloc,           // tp_dealloc
  0,                         // tp_print
  0,                         // tp_getattr
  0,                         // tp_setattr
  0,                         // tp_compare
  0,                         // tp_repr
  0,                         // tp_as_number
  0,                         // tp_as_sequence
  0,                         // tp_as_mapping
  0,                         // tp_hash
  0,                         // tp_call
  0,                         // tp_str
  0,                         // tp_getattro
  0,                         // tp_setattro
  0,                         // tp_as_buffer
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, // tp_flags
  "Filter objects",            // tp_doc
  0,                         // tp_traverse
  0,                         // tp_clear
  0,                         // tp_richcompare
  0,                         // tp_weaklistoffset
  0,                         // tp_iter
  0,                         // tp_iternext
  PyFilter_methods,                         // tp_methods
  0,                         // tp_members
  0,                         // tp_getset
  0,                         // tp_base
  0,                         // tp_dict
  0,                         // tp_descr_get
  0,                         // tp_descr_set
  0,                         // tp_dictoffset
  PyFilter::init,              // tp_init
  0,                         // tp_alloc
  PyFilter::make_new,        // tp_new
};

int PyFilter::init(PyObject* self, PyObject* args, PyObject* kwds) {
  Real samplerate;
  int channels;

  PyObject *inputA;
  PyObject *inputB;

  PyArrayObject *coeffsA;
  PyArrayObject *coeffsB;

  // Parse the arguments
  if (!PyArg_ParseTuple(args, "OOfi", &inputA, &inputB, &samplerate, &channels)){
    PyErr_SetString(PyExc_ValueError,
                    "must pass two arrays of coefficients (a and b) and the samplerate and channels as arguments");
    return -1;
  }

  coeffsA = (PyArrayObject *) PyArray_ContiguousFromObject(inputA, PyArray_FLOAT, 2, 2);
  coeffsB = (PyArrayObject *) PyArray_ContiguousFromObject(inputB, PyArray_FLOAT, 2, 2);

  if ((coeffsA == NULL) || (coeffsB == NULL)) {
    PyErr_SetString(PyExc_ValueError,
                    "the coefficients arrays must be of type float");
    
    return -1;
  }

  // Check that the coefficients have same number of columns as the channels passed as parameter
  if ((coeffsA->dimensions[1] != channels) || (coeffsB->dimensions[1] != channels)) {
    PyErr_SetString(PyExc_ValueError,
                    "the coefficients arrays must have as many columns as channels");
    
    return -1;
  }

  int a_rows = coeffsA->dimensions[0];
  int a_cols = coeffsA->dimensions[1];
  Real* a_data = (Real*)PyArray_DATA(coeffsA);

  MatrixXR a = Eigen::Map<MatrixXRscipy>(a_data, a_rows, a_cols);

  int b_rows = coeffsB->dimensions[0];
  int b_cols = coeffsB->dimensions[1];
  Real* b_data = (Real*)PyArray_DATA(coeffsB);

  MatrixXR b = Eigen::Map<MatrixXRscipy>(b_data, b_rows, b_cols);


  if(((PyFilter*)self)->base != NULL)
    delete ((PyFilter*)self)->base;
  
  // Create the filter and set it up
  ((PyFilter*)self)->base = new Filter::Filter(b, a, samplerate, channels);
  ((PyFilter*)self)->base->setup();
  
  return 0;
}

PyObject* PyFilter::process(PyObject* self, PyObject* args){
  PyObject *input;
  PyArrayObject *in_array;

  // check that the input argument is an array
  if (!PyArg_ParseTuple(args, "O", &input)){
    PyErr_SetString(PyExc_ValueError,
                    "process only takes one argument of type array");

    return NULL;
  }

  in_array = (PyArrayObject *) PyArray_ContiguousFromObject(input, PyArray_FLOAT, 2, 2);
  
  if (in_array == NULL){
    PyErr_SetString(PyExc_ValueError,
                    "array must be of type float");
    
    return NULL;
  }


  // check that the input array has 2 dimensions
  if (in_array->nd != 2 || in_array->descr->type_num != PyArray_FLOAT) { 
    PyErr_SetString(PyExc_ValueError,
                    "array must be two-dimensional and of type float");
    
    return NULL;
  }

  // prepare the input array
  int in_samples = in_array->dimensions[0];
  int in_channels = in_array->dimensions[1];

  // prepare the dimensions variables
  int out_samples = in_samples;
  int out_channels = ((PyFilter*)self)->base->channels();
  
  // check that the input array has the right number of channels
  if (!((in_channels == 1) || (in_channels == out_channels))) {
    PyErr_SetString(PyExc_ValueError,
                    "the input must be an array of one single channel (1 column) or as many channels as filters");
    
    return NULL;

  }
  
  // prepare resulting array
  int dims[] = {out_samples, out_channels};
  PyObject * out_array = PyArray_FromDims(2, dims, PyArray_FLOAT);

  if (out_array == NULL)
    return NULL;

  // call the process method
  int in_rows = in_array->dimensions[0];
  int in_cols = in_array->dimensions[1];
  Real* in_data = (Real*)PyArray_DATA(in_array);

  MatrixXR in_matrix = Eigen::Map<MatrixXRscipy>(in_data, in_rows, in_cols);
  MatrixXR out_matrix(in_rows, out_channels);

  ((PyFilter*)self)->base->process(in_matrix, &out_matrix);


  int out_rows = ((PyArrayObject*)out_array)->dimensions[0];
  int out_cols = ((PyArrayObject*)out_array)->dimensions[1];
  Real* out_data = (Real*)PyArray_DATA(out_array);

  Eigen::Map<MatrixXRscipy>(out_data, out_rows, out_cols) = out_matrix;
  
  return out_array;
}
