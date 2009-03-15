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

#ifndef CORRELATION_H
#define CORRELATION_H

#include "Typedefs.h"
#include "Debug.h"

#include "FFT.h"
#include "IFFT.h"

#include <limits>

class Correlation {
protected:
  // Internal parameters
  int _inputLengthA;
  int _inputLengthB;
  int _minLag;
  int _maxLag;
  bool _useFFT;
  int _fftSize;

  // Internal variables
  FFT _fft;
  IFFT _ifft;

  MatrixXC _fftA;
  MatrixXC _fftB;
  MatrixXR _result;

public:
  Correlation(int inputLengthA, int inputLengthB, Real maxLag = std::numeric_limits<Real>::infinity(), Real minLag = -std::numeric_limits<Real>::infinity());
  Correlation(int inputLengthA, int inputLengthB, Real maxLag, Real minLag, bool useFFT);

  ~Correlation();

  void setup();

  void process(const MatrixXR& inputA, const MatrixXR& inputB, MatrixXR* autocorrelation);

  void reset();
};

#endif  /* CORRELATION_H */
