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

#include "pymfcc.h"

#define NO_IMPORT_ARRAY
#include "ndarrayobject.h"

#include "typedefs.h"

using namespace std;

PyMethodDef PyMFCC_methods[] = {
  { "process",            (PyCFunction)PyMFCC::process, METH_VARARGS,
    "MFCC.process(input_array) processes the input samples \"input_array\"" },
  { NULL }  /* Sentinel */
};

PyTypeObject PyMFCCType = {
  PyObject_HEAD_INIT(NULL)
  0,                         // ob_size
  "cricaudio.MFCC",           // tp_name
  sizeof(PyMFCC),            // tp_basicsize
  0,                         // tp_itemsize
  PyMFCC::dealloc,           // tp_dealloc
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
  "MFCC objects",            // tp_doc
  0,                         // tp_traverse
  0,                         // tp_clear
  0,                         // tp_richcompare
  0,                         // tp_weaklistoffset
  0,                         // tp_iter
  0,                         // tp_iternext
  PyMFCC_methods,                         // tp_methods
  0,                         // tp_members
  0,                         // tp_getset
  0,                         // tp_base
  0,                         // tp_dict
  0,                         // tp_descr_get
  0,                         // tp_descr_set
  0,                         // tp_dictoffset
  PyMFCC::init,              // tp_init
  0,                         // tp_alloc
  PyMFCC::make_new,          // tp_new
};

int PyMFCC::init(PyObject* self, PyObject* args, PyObject* kwds) {
  Real lowFreq;
  Real highFreq;
  int numBands;
  Real samplerate;
  int spectrumLength;
  int numCoeffs;
  Real minSpectrum;
  Real power;


  // default constructor with no argument
  if (!PyArg_ParseTuple(args, "ffifiiff", &lowFreq, &highFreq, &numBands, &samplerate, &spectrumLength, &numCoeffs, &minSpectrum, &power)){
    PyErr_SetString(PyExc_ValueError,
                    "the arguments must be (float)lowFreq, (float)highFreq, (int)numBands, (float)samplerate, (int)spectrumLength, (int)numCoeffs, (float)minSpectrum, (float)power");
    
    return -1;
  }

  if(((PyMFCC*)self)->base != NULL)
    delete ((PyMFCC*)self)->base;
      
  ((PyMFCC*)self)->base = new MFCC::MFCC(lowFreq, highFreq, numBands, samplerate, spectrumLength, numCoeffs, minSpectrum, power);
  ((PyMFCC*)self)->base->setup();

  return 0;
}

PyObject* PyMFCC::process(PyObject* self, PyObject* args){
  PyObject *input;
  PyArrayObject *in_array;

  // check that the input argument is an array
  if (!PyArg_ParseTuple(args, "O", &input))
    return NULL;

  in_array = (PyArrayObject *) PyArray_ContiguousFromObject(input, PyArray_FLOAT, 2, 2);
  
  if (in_array == NULL)
    return NULL;


  // check that the input array has 2 dimensions
  if (in_array->nd != 2 || in_array->descr->type_num != PyArray_FLOAT) { 
    PyErr_SetString(PyExc_ValueError,
                    "array must be two-dimensional and of type float");
    
    return NULL;
  }

  // prepare the input array
  int in_rows = in_array->dimensions[0];
  int in_cols = in_array->dimensions[1];
  Real* in_data = (Real*)PyArray_DATA(in_array);

  MatrixXR in_matrix = Eigen::Map<MatrixXRscipy>(in_data, in_rows, in_cols);
  
  // TODO: move this check in the mfcc.cpp
  // check that the input array has the right number of channels
  if (in_cols != 1) {
    PyErr_SetString(PyExc_ValueError,
                    "the number of columns of the input array must be 1");
    
    return NULL;

  }

  // prepare resulting array
  int numCoeffs = ((PyMFCC*)self)->base->numCoeffs();
  int dims[] = {numCoeffs, 1};
  PyObject* out_array = PyArray_FromDims(2, dims, PyArray_FLOAT);
  
  if (out_array == NULL)
    return NULL;

  MatrixXR out_matrix(numCoeffs, 1);

  ((PyMFCC*)self)->base->process(in_matrix, &out_matrix);

  int out_rows = ((PyArrayObject*)out_array)->dimensions[0];
  int out_cols = ((PyArrayObject*)out_array)->dimensions[1];
  Real* out_data = (Real*)PyArray_DATA(out_array);

  Eigen::Map<MatrixXRscipy>(out_data, out_rows, out_cols) = out_matrix;
  
  return out_array;
}
