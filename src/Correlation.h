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

/**
  * @class Correlation
  *
  * @brief Correlation algorithm.
  *
  * This class represents an object to perform a correlation of two vectors.
  *
  * The correlation can be performed using two methods:
  * -# Direct method 
  * -# FFT method
  *
  * The Direct method performs faster than the FFT method on vectors of small sizes.
  * The critical sizes of the vectors depends on the platform.
  *
  * @author Ricard Marxer
  *
  * @sa Autocorrelation
  */
class Correlation {
protected:
  // Internal parameters
  int _inputSizeA;
  int _inputSizeB;
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
  /**
     Constructs an Autocorrelation object with the specified @a inputSize, 
     @a maxLag and @a minLag settings.
     
     @param inputSize size of the inputs arrays to be autocorrelated,
     must be > 0.
     The algorithm performs faster for sizes which are a power of 2.
     
     @param maxLag maximum lag to be calculated
     @param minLag minimum lag to be calculated
     @param useFFT determines whether or not to use the FFT method
  */
  Correlation(int inputSizeA, int inputSizeB, int maxLag = std::numeric_limits<int>::max(), int minLag = -std::numeric_limits<int>::max());
  Correlation(int inputSizeA, int inputSizeB, int maxLag, int minLag, bool useFFT);

  /**
     Destroys the Correlation algorithm and frees its resources.
  */
  ~Correlation();

  void setup();
  void reset();

  /**
     Performs a Correlation between each of the rows of @a framesA and
     each of the rows of @b framesB respectively.
     Puts the resulting correlations in the rows of @a correlation.
     
     @param framesA matrix of Real values.  The number of columns of @a 
     framesA must be equal to the inputSizeA property.

     @param framesB matrix of Real values.  The number of columns of @a 
     framesB must be equal to the inputSizeB property.

     Note that @a framesA and @a framesB should have the same number of rows.

     @param correlation pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows as @a framesA (and @a framesB) and maxLag - minLag columns.

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& framesA, const MatrixXR& framesB, MatrixXR* autocorrelation);

  /**
     Returns the size of the first input arrays to be correlated.
     The default is 1024.
     
     @sa inputSizeB(), setInputSizeA(), setInputSizeB()
  */
  int inputSizeA() const;  

  /**
     Specifies the @a size of the first of the input arrays.
     The given @a size must be higher than 0.
     Note that if @a size is a power of 2 the algorithm will perform faster.
     
     @sa inputSizeA(), inputSizeB(), setInputSizeB()
  */
  void setInputSizeA( int size, bool callSetup = true );

  /**
     Returns the size of the second of the input arrays to be correlated.
     The default is 1024.
     
     @sa inputSizeA(), setInputSizeA(), setInputSizeB()
  */
  int inputSizeB() const;  

  /**
     Specifies the @a size of the second of the input arrays.
     The given @a size must be higher than 0.
     Note that if @a size is a power of 2 the algorithm will perform faster.
     
     @sa inputSizeA(), inputSizeB(), setInputSizeA()
  */
  void setInputSizeB( int size, bool callSetup = true );

  /**
     Returns the minimum lag to be calculated.
     The default is -max(_inputSizeA, _inputSizeB) + 1.
     
     @sa maxLag(), setMinLag(), setMaxLag()
  */
  int minLag() const;  
  
  /**
     Specifies the minimum @a lag of the ocorrelation.
     The given @a lag will be constratined between - max( inputSizeA, inputSizeB ) + 1 and min( inputSizeA, inputSizeB ).
     Note that the lag should be smaller than the maxLag.

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
     Specifies the maximum @a lag of the correlation.
     The given @a lag will be constratined between - max( inputSizeA, inputSizeB ) + 1 and min( inputSizeA, inputSizeB ).

     Note that the lag should be larger than the maxLag.

     @sa minLag(), maxLag(), setMinLag()
  */
  void setMaxLag( int lag, bool callSetup = true );

  /**
     Returns @c true if the FFT method should be used for the correlation.
     The default is True for inputSize larger than 128; otherwise it is False.
     
     @sa setUseFFT()
  */
  bool useFFT() const;  
  
  /**
     Specifies whether the autocorrelation should be performed using the FFT method.
     
     @sa useFFT()
  */
  void setUseFFT( bool useFFT, bool callSetup = true );
};

#endif  /* CORRELATION_H */
