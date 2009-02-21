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

#ifndef AUTOCORRELATION_H
#define AUTOCORRELATION_H

#include "Typedefs.h"
#include "Debug.h"

#include "FFT.h"
#include "IFFT.h"

#include <limits>

class Autocorrelation {
protected:
  // Internal parameters
  int _inputLength;
  int _minLag;
  int _maxLag;
  bool _useFFT;

  // Internal variables
  FFT _fft;
  IFFT _ifft;

public:
  Autocorrelation(int inputLength, int maxLag = std::numeric_limits<Real>::infinity(), int minLag = 0);
  Autocorrelation(int inputLength, int maxLag, int minLag, bool useFFT);

  ~Autocorrelation();

  void setup();

  void process(const MatrixXR& input, MatrixXR* autocorrelation);

  void reset();
};

#endif  /* AUTOCORRELATION_H */