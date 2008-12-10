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

static PyObject* CRicaudio__Module = NULL;

static PyMethodDef CRicaudio__Methods[] = {
  { NULL,           NULL,                     0,           NULL } // Sentinel
};


extern PyTypeObject PyMeddisType;
extern PyTypeObject PyMFCCType;
extern PyTypeObject PyFilterType;
extern PyTypeObject PyAOKType;

extern "C" void initricaudio() {
  
  // import our wrapper types
  if ((PyType_Ready(&PyMeddisType) < 0) || (PyType_Ready(&PyMFCCType) < 0) || (PyType_Ready(&PyFilterType) < 0) || (PyType_Ready(&PyAOKType) < 0)) {
    
    cerr << "Unable to instantiate CRicaudio's wrapper types." << endl;
    return;
  }

  CRicaudio__Module = Py_InitModule3("ricaudio", CRicaudio__Methods,
                                     "Module that allows access to cricaudio plugins and algorithms.");
    
  if (CRicaudio__Module == NULL)
    return;
  
  // import the NumPy C api
  int numpy_error = _import_array();
  if (numpy_error) {
    cerr << "Unable to import NumPy C API from CRicaudio module. Error code = " << numpy_error << endl;
    return;
  }

  // insert the classes
  Py_INCREF(&PyMeddisType);
  PyModule_AddObject(CRicaudio__Module, (char*)"Meddis", (PyObject*)&PyMeddisType);

  Py_INCREF(&PyMFCCType);
  PyModule_AddObject(CRicaudio__Module, (char*)"MFCC", (PyObject*)&PyMFCCType);

  Py_INCREF(&PyFilterType);
  PyModule_AddObject(CRicaudio__Module, (char*)"Filter", (PyObject*)&PyFilterType);

  Py_INCREF(&PyAOKType);
  PyModule_AddObject(CRicaudio__Module, (char*)"AOK", (PyObject*)&PyAOKType);

}
