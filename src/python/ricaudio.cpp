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

#include <Python.h>
#include <iostream>

using namespace std;

// numpy import
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_API
#include "arrayobject.h"

#include "pymeddis.h"
#include "pymfcc.h"
#include "pyfilter.h"
#include "pyaok.h"
#include "pyfft.h"
#include "pywindow.h"
#include "pyspectralreassignment.h"

static PyObject* Ricaudio__Module = NULL;

static PyMethodDef Ricaudio__Methods[] = {
  { NULL,           NULL,                     0,           NULL } // Sentinel
};


extern PyTypeObject PyMeddisType;
extern PyTypeObject PyMFCCType;
extern PyTypeObject PyFilterType;
extern PyTypeObject PyAOKType;
extern PyTypeObject PyFFTType;
extern PyTypeObject PyWindowType;
extern PyTypeObject PySpectralReassignmentType;

PyMODINIT_FUNC
initricaudio(void) {
  
  // import our wrapper types
  if ((PyType_Ready(&PyMeddisType) < 0) || \
      (PyType_Ready(&PyMFCCType) < 0) || \
      (PyType_Ready(&PyFilterType) < 0) || \
      (PyType_Ready(&PyAOKType) < 0) || \
      (PyType_Ready(&PyFFTType) < 0) || \
      (PyType_Ready(&PySpectralReassignmentType) < 0) || \
      (PyType_Ready(&PyWindowType) < 0)) {
    
    cerr << "Unable to instantiate Ricaudio's wrapper types." << endl;
    return;
  }

  Ricaudio__Module = Py_InitModule3("ricaudio", Ricaudio__Methods,
                                    "Module that allows access to ricaudio plugins and algorithms.");
    
  if (Ricaudio__Module == NULL)
    return;
  
  // import the NumPy C api
  int numpy_error = _import_array();
  if (numpy_error) {
    cerr << "Unable to import NumPy C API from Ricaudio module. Error code = " << numpy_error << endl;
    return;
  }

  // insert the classes
  Py_INCREF(&PyMeddisType);
  PyModule_AddObject(Ricaudio__Module, (char*)"Meddis", (PyObject*)&PyMeddisType);

  Py_INCREF(&PyMFCCType);
  PyModule_AddObject(Ricaudio__Module, (char*)"MFCC", (PyObject*)&PyMFCCType);

  Py_INCREF(&PyFilterType);
  PyModule_AddObject(Ricaudio__Module, (char*)"Filter", (PyObject*)&PyFilterType);

  Py_INCREF(&PyAOKType);
  PyModule_AddObject(Ricaudio__Module, (char*)"AOK", (PyObject*)&PyAOKType);

  Py_INCREF(&PyFFTType);
  PyModule_AddObject(Ricaudio__Module, (char*)"FFT", (PyObject*)&PyFFTType);

  Py_INCREF(&PyWindowType);
  PyModule_AddObject(Ricaudio__Module, (char*)"Window", (PyObject*)&PyWindowType);
  
  Py_INCREF(&PySpectralReassignmentType);
  PyModule_AddObject(Ricaudio__Module, (char*)"SpectralReassignment", (PyObject*)&PySpectralReassignmentType);

}
