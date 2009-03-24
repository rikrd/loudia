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
  * @brief Algorithm to perform an autocorrelation of vectors of Real values.
  *
  * This class represents an object to perform a correlation of a vector with itself.
  *
  * The correlation can be performed using two methods:
  * -# Direct method 
  * -# FFT method
  *
  * The Direct method consists in applying the correlation formula directly
  * in the time domain.
  *
  * The FFT method consists in performing an Fast Fourier Transform of the
  * vector and multiply it by its conjugate.
  * Finally the algorithm applies an IFFT to the result of the 
  * multiplication in order to obtain the autocorrelation for all
  * the time lags.
  *
  * The Direct method performs faster than the FFT method only
  * on vectors of small sizes. The decision point for selecting one of
  * the two methods depends on the platform.
  *
  * The method performed can be specified using setUseFft().
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
  bool _useFft;

  // Internal variables
  int _calcMinLag;
  int _calcMaxLag;
  FFT _fft;
  IFFT _ifft;
  MatrixXC _tempFft;
  MatrixXR _temp;

public:
  /**
     Constructs an Autocorrelation object with the specified @a inputSize, 
     @a maxLag and @a minLag settings.
     
     @param inputSize size of the inputs arrays to be autocorrelated,
     must be > 0.
     The algorithm performs faster for sizes which are a power of 2.
     
     @param maxLag maximum lag to be calculated
     @param minLag minimum lag to be calculated
     @param useFft determines whether or not to use the FFT method
  */
  Autocorrelation(int inputSize, int maxLag, int minLag, bool useFft);
  Autocorrelation(int inputSize, int maxLag, int minLag);
  Autocorrelation(int inputSize, int maxLag);
  Autocorrelation(int inputSize = 1024);

  /**
     Destroys the Autocorrelation algorithm and frees its resources.
  */
  ~Autocorrelation();

  void setup();
  void reset();

  /**
     Performs an autocorrelation on each of the rows of @a frames.
     Puts the resulting autocorrelations in the rows of @a autocorrelation.
     
     @param frames matrix of Real values.  The number of columns of @a 
     frames must be equal to the inputSize property.
     
     @param autocorrelation pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows as @a frames and maxLag - minLag columns.

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& frames, MatrixXR* autocorrelation);

  /**
     Returns the size of the input arrays to be autocorrelated.
     The default is 1024.
     
     @sa setInputSize()
  */
  int inputSize() const;  

  /**
     Specifies the @a size of the input.
     The given @a size must be higher than 0.
     Note that if @a size is a power of 2 the algorithm will perform faster.
     
     @sa inputSize()
  */
  void setInputSize( int size, bool callSetup = true );

  /**
     Returns the minimum lag to be calculated.
     The default is 0.
     
     @sa maxLag(), setMinLag(), setMaxLag()
  */
  int minLag() const;  
  
  /**
     Specifies the minimum @a lag of the autocorrelation.
     The given @a lag will be constratined between -inputSize + 1 and inputSize.
     
     @sa minLag(), maxLag(), setMaxLag()
  */  
  void setMinLag( int lag, bool callSetup = true );

  /**
     Returns the maximum lag to be calculated.
     The default is inputSize.
     
     @sa minLag(), setMinLag(), setMaxLag()
  */
  int maxLag() const;  
  
  /**
     Specifies the maximum @a lag of the autocorrelation.
     The given @a lag will be constratined between -inputSize + 1 and inputSize.
     
     @sa minLag(), maxLag(), setMinLag()
  */
  void setMaxLag( int lag, bool callSetup = true );

  /**
     Returns @c true if the FFT method should be used for the autocorrelation.
     The default is True for inputSize larger than 128; otherwise it is False.
     
     @sa setUseFft()
  */
  bool useFft() const;  
  
  /**
     Specifies whether the autocorrelation should be performed using the FFT method.
     
     @sa useFft()
  */
  void setUseFft( bool useFft, bool callSetup = true );
};

#endif  /* AUTOCORRELATION_H */
