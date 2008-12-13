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

#include "Eigen/Core"
#include "Eigen/Array"

#include "pyfft.h"

#define NO_IMPORT_ARRAY
#include "ndarrayobject.h"

#include "typedefs.h"

using namespace std;

PyMethodDef PyFFT_methods[] = {
  { "process",            (PyCFunction)PyFFT::process, METH_VARARGS,
    "FFT.process(input_array) processes the input samples \"input_array\"" },
  { NULL }  /* Sentinel */
};

PyTypeObject PyFFTType = {
  PyObject_HEAD_INIT(NULL)
  0,                         // ob_size
  "cricaudio.FFT",           // tp_name
  sizeof(PyFFT),            // tp_basicsize
  0,                         // tp_itemsize
  PyFFT::dealloc,           // tp_dealloc
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
  "FFT objects",            // tp_doc
  0,                         // tp_traverse
  0,                         // tp_clear
  0,                         // tp_richcompare
  0,                         // tp_weaklistoffset
  0,                         // tp_iter
  0,                         // tp_iternext
  PyFFT_methods,                         // tp_methods
  0,                         // tp_members
  0,                         // tp_getset
  0,                         // tp_base
  0,                         // tp_dict
  0,                         // tp_descr_get
  0,                         // tp_descr_set
  0,                         // tp_dictoffset
  PyFFT::init,              // tp_init
  0,                         // tp_alloc
  PyFFT::make_new,          // tp_new
};

int PyFFT::init(PyObject* self, PyObject* args, PyObject* kwds) {
  int frameSize;
  int fftSize;
  bool zeroPhase;

  // default constructor with no argument
  if (!PyArg_ParseTuple(args, "iib", &frameSize, &fftSize, &zeroPhase)){
    PyErr_SetString(PyExc_ValueError,
                    "the arguments must be (int)frameSize, (int)fftSize, (bool)zeroPhase");
    
    return -1;
  }

  if(((PyFFT*)self)->base != NULL)
    delete ((PyFFT*)self)->base;
      
  ((PyFFT*)self)->base = new FFT::FFT(frameSize, fftSize, zeroPhase);
  ((PyFFT*)self)->base->setup();

  return 0;
}

PyObject* PyFFT::process(PyObject* self, PyObject* args){
  PyObject *input;
  PyArrayObject *in_array;

  // check that the input argument is an array
  if (!PyArg_ParseTuple(args, "O", &input))
    return NULL;

  in_array = (PyArrayObject *) PyArray_ContiguousFromObject(input, PyArray_CFLOAT, 2, 2);
  
  if (in_array == NULL)
    return NULL;


  // check that the input array has 2 dimensions
  if (in_array->nd != 2 || in_array->descr->type_num != PyArray_CFLOAT) { 
    PyErr_SetString(PyExc_ValueError,
                    "array must be two-dimensional and of type float");
    
    return NULL;
  }

  // prepare the input array
  int in_rows = in_array->dimensions[0];
  int in_cols = in_array->dimensions[1];
  Complex* in_data = (Complex*)PyArray_DATA(in_array);

  MatrixXC in_matrix = Eigen::Map<MatrixXCscipy>(in_data, in_rows, in_cols);
  
  // prepare resulting array
  int fftSize = ((PyFFT*)self)->base->fftSize();
  int dims[] = {in_rows, fftSize};
  PyObject* out_array = PyArray_FromDims(2, dims, PyArray_CFLOAT);
  
  if (out_array == NULL)
    return NULL;

  MatrixXC out_matrix(in_rows, fftSize);

  // process the data
  ((PyFFT*)self)->base->process(in_matrix, &out_matrix);

  int out_rows = ((PyArrayObject*)out_array)->dimensions[0];
  int out_cols = ((PyArrayObject*)out_array)->dimensions[1];
  Complex* out_data = (Complex*)PyArray_DATA(out_array);

  // put the data in the resulting array
  Eigen::Map<MatrixXCscipy>(out_data, out_rows, out_cols) = out_matrix;
  
  return out_array;
}
