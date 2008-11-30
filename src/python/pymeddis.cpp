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

#include "Eigen/Core"
#include "Eigen/Array"

#include "pymeddis.h"

#define NO_IMPORT_ARRAY
#include "ndarrayobject.h"

#include "typedefs.h"

using namespace std;

PyMethodDef PyMeddis_methods[] = {
  { "process",            (PyCFunction)PyMeddis::process, METH_VARARGS,
    "Meddis.process(input_array) processes the input samples \"input_array\"" },
  { NULL }  /* Sentinel */
};

PyTypeObject PyMeddisType = {
  PyObject_HEAD_INIT(NULL)
  0,                         // ob_size
  "cricaudio.Meddis",           // tp_name
  sizeof(PyMeddis),            // tp_basicsize
  0,                         // tp_itemsize
  PyMeddis::dealloc,           // tp_dealloc
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
  "Meddis objects",            // tp_doc
  0,                         // tp_traverse
  0,                         // tp_clear
  0,                         // tp_richcompare
  0,                         // tp_weaklistoffset
  0,                         // tp_iter
  0,                         // tp_iternext
  PyMeddis_methods,                         // tp_methods
  0,                         // tp_members
  0,                         // tp_getset
  0,                         // tp_base
  0,                         // tp_dict
  0,                         // tp_descr_get
  0,                         // tp_descr_set
  0,                         // tp_dictoffset
  PyMeddis::init,              // tp_init
  0,                         // tp_alloc
  PyMeddis::make_new,          // tp_new
};

int PyMeddis::init(PyObject* self, PyObject* args, PyObject* kwds) {
  Real samplerate;
  int channels;
  
  // default constructor with no argument
  if (!PyArg_ParseTuple(args, "fi", &samplerate, &channels)){
    PyErr_SetString(PyExc_ValueError,
                    "the arguments must be a float and an integer");
    
    return -1;
  }

  if(((PyMeddis*)self)->base != NULL)
    delete ((PyMeddis*)self)->base;
      
  ((PyMeddis*)self)->base = new Meddis::Meddis(samplerate, channels);
  ((PyMeddis*)self)->base->setup();

  return 0;
}

PyObject* PyMeddis::process(PyObject* self, PyObject* args){
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
  
  // TODO: move this check in the meddis.cpp
  // check that the input array has the right number of channels
  if (in_cols != ((PyMeddis*)self)->base->channels()) {
    PyErr_SetString(PyExc_ValueError,
                    "the number of channels must be the same as the one set for the Meddis object");
    
    return NULL;

  }

  // prepare resulting array
  int dims[] = {in_rows, in_cols};
  PyObject* out_array = PyArray_FromDims(2, dims, PyArray_FLOAT);
  
  if (out_array == NULL)
    return NULL;

  MatrixXR out_matrix(in_rows, in_cols);

  ((PyMeddis*)self)->base->process(in_matrix, &out_matrix);

  int out_rows = ((PyArrayObject*)out_array)->dimensions[0];
  int out_cols = ((PyArrayObject*)out_array)->dimensions[1];
  Real* out_data = (Real*)PyArray_DATA(out_array);

  Eigen::Map<MatrixXRscipy>(out_data, out_rows, out_cols) = out_matrix;
  
  //return PyArray_Return(out_array);
  return out_array;
}
