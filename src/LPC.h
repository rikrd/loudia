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

#ifndef LPC_H
#define LPC_H

#include "Typedefs.h"
#include "Debug.h"

#include "Filter.h"
#include "Autocorrelation.h"

/**
  * @class LPC
  *
  * @brief Algorithm to calculate the Linear Predictive Coding of vectors of Real values.
  *
  * This class represents an object to perform a Linear Predictive Coding on vectors 
  * of Real values.  Which is a useful technique for estimating a parametric representation
  * of a spectrum magnitude.  The algorithm estimates a set of M coefficients of a IIR filter
  * whose frequency response approximates the vector of Reals passed as input.
  *
  * This algorithm implements the Levinson-Durbin recursion for solving the
  * following linear equation system:
  *
  * R a = r
  *
  * where R is the Topelitz matrix made of the first M - 1 autocorrelation 
  * coefficients of the input vector and r is a vector made of M - 1 
  * autocorrelation coefficients starting from the second of the input
  * vector.
  *
  * Optionally a pre-emphasis FIR filter may be applied to the input vector
  * in order to enhance estimation of higher frequencies. The pre-emphasis filter
  * consists of a 2 coefficient filter of the form b = [1, b1] where usually:
  *
  * 0.96 <= b1 <= 0.99
  *
  * The b1 coefficient defaults to 0, but can be specified using setPreEmphasis().
  *
  * @author Ricard Marxer
  *
  * @sa MelBands, Bands, MFCC
  */
class LPC {
protected:
  // Internal parameters
  int _inputSize;
  int _coefficientCount;
  Real _preEmphasis;
  
  // Internal variables
  MatrixXR _pre;
  MatrixXR _preRow;
  MatrixXR _temp;
  MatrixXR _acorr;

  Filter _preFilter;
  Autocorrelation _acorrelation;

public:
  /**
     Constructs an LPC object with the specified @a inputSize, 
     @a coefficientCount and @a preEmphasis settings.
     
     @param inputSize size of the inputs arrays,
     must be > 0.
     The algorithm performs faster for sizes which are a power of 2.
     
     @param coefficientCount number of coefficients to be estimated
     @param preEmphasis second coefficient of the FIR pre-emphasis filter
  */
  LPC(int inputSize = 1024, int coefficientCount = 15, Real preEmphasis = 0.0);

  /**
     Destroys the algorithm and frees its resources.
  */
  ~LPC();

  void setup();
  void reset();

  /**
     Performs an LPC on each of the rows of @a frames.
     Puts the resulting LPC coefficients in the rows of @a lpcCoefficients, the
     reflection coefficients in @a reflectionCoefficients and the error in @a error.
     
     @param frames matrix of Real values.  The number of columns of @a 
     frames must be equal to the input size specified using setInputSize().
     
     @param lpcCoefficients pointer to a matrix of Real values for the LPC coefficients.
     The matrix should have the same number of rows as @a frames and coefficientCount columns.

     @param reflectionCoefficients pointer to a matrix of 
     Real values for the reflection coefficients.
     The matrix should have the same number of rows as @a frames 
     and coefficientCount + 1 columns.
     
     @param error pointer to a matrix of 
     Real values for the LPC error gain.
     The matrix should have the same number of rows as @a frames 
     and 1 single column.

     Note that if the output matrices are not of the required sizes they will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& frames, MatrixXR* lpcCoefficients,
               MatrixXR* reflectionCoefficients, MatrixXR* error);

  /**
     Returns the size of the input arrays.
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
     Returns the number of coefficients to be calculated.
     The default is 15.
     
     @sa setCoefficientCount()
  */
  int coefficientCount() const;

  /**
     Specifies the @a count of coefficients to be calculated.
     The given @a count must be in the range between 0 and (input size - 1).
          
     @sa coefficientCount()
  */
  void setCoefficientCount( int count, bool callSetup = true );

  /**
     Returns the second coefficient of the FIR preemphasis filter.
     The default is 0.0.
     
     @sa setPreEmphasis()
  */
  Real preEmphasis() const;

  /**
     Specifies the second @a coefficient of the FIR preemphasis filter.
     
     @sa preEmphasis()
  */
  void setPreEmphasis( Real coefficient, bool callSetup = true );

};

#endif  /* LPC_H */
