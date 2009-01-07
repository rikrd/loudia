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

#ifndef ODFCOMPLEX_H
#define ODFCOMPLEX_H

#include "typedefs.h"
#include "debug.h"

#include "window.h"
#include "fft.h"
#include "unwrap.h"

class ODFComplex {
protected:
  // Internal parameters
  int _frameLength;
  int _fftLength;
  bool _zeroPhase;
  Window::WindowType _windowType;
  
  // Internal variables
  Window _window;
  FFT _fft;
  Unwrap _unwrap;

  MatrixXR _windowed;
  MatrixXC _ffted;
  MatrixXC _spectrum;
  MatrixXR _unwrappedAngle;

public:
  ODFComplex(int frameLength, int fftLength, Window::WindowType windowType = Window::RECTANGULAR, bool zeroPhase = true);

  ~ODFComplex();

  void setup();

  void process(MatrixXR samples, MatrixXR* odfValue);

  void reset();

};

#endif  /* ODFCOMPLEX_H */
