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

/**
  * @class Autocorrelation
  *
  * @brief Autocorrelation algorithm.
  *
  * This class represents an object to perform a correlation of a vector with itself.
  *
  * The correlation can be performed using two methods:
  * -# Direct method 
  * -# FFT method
  *
  * The Direct method performs faster than the FFT method on vectors of small sizes.
  * The critical size of the vector depends on the platform.
  *
  * @author Ricard Marxer
  *
  * @sa Correlation, PitchACF
  */
class Autocorrelation {
protected:
  // Internal parameters
  int _inputSize;
  int _minLag;
  int _maxLag;
  bool _useFFT;

  // Internal variables
  FFT _fft;
  IFFT _ifft;

public:
  Autocorrelation(int inputSize = 1024, int maxLag = std::numeric_limits<int>::max(), int minLag = 0);
  Autocorrelation(int inputSize, int maxLag, int minLag, bool useFFT);

  ~Autocorrelation();

  void setup();
  void reset();

  void process(const MatrixXR& input, MatrixXR* autocorrelation);

  int inputSize() const;  
  void setInputSize( int size );

  int minLag() const;  
  void setMinLag( int lag );

  int maxLag() const;  
  void setMaxLag( int lag );

  bool useFFT() const;  
  void setUseFFT( bool useFFT );
};

#endif  /* AUTOCORRELATION_H */
