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

#ifndef FILTER_H
#define FILTER_H

#include "Typedefs.h"
#include "Debug.h"

/**
  * @class Filter
  *
  * @brief Algorithm to apply one or several IIR filters given the Real value coefficients.
  *
  * This class represents an object to apply one or several IIR filters.
  * The coefficients must be manually set by the user.  To create and use
  * some special parametrized filters such as Low Pass, High Pass, Band Pass and Band Stop
  * filters refer to @l IIRFilter.
  *
  * For Complex coefficient filters, refer to @l FilterComplex.
  *
  * This filter implementation allows single and multiple channel filtering.
  * The number of channels is defined by the @l channels property.
  * The a and b coefficients of the filter are specified by two matrices of
  * Real values. The rows of the matrix are the time indices of the filter
  * and the columns (if more than one) are the channels.
  *
  * Three different situations are possible, depending on the number of columns
  * in the coefficients matrices and in the input matrix:
  * - if the number of columns of the coefficients matrices 
  * are equal to the number columns in the input matrix, then
  * each column of the input matrix is filtered by a column 
  * of the coefficients matrix. This is the situation when trying to
  * filter differently all the channels of a multi-channel signal.
  * - if the coefficients matrix has one single column 
  * and the input matrix has multiple columns then each 
  * column of the input matrix is filtered by the single column of the
  * coefficients matrices.  This is the situation when trying to
  * filter equally all the channels of a multi-channel signal.
  * - if the coefficients matrices have multiple columns each and the 
  * input matrix has multiple columns then the column of the input matrix is filtered 
  * by each of the columns of the coefficients matrices.  This is the situation
  * when applying a filterbank to a single channel signal.
  *
  * Note that in all cases the number of columns in a and b coefficients matrices
  * must be the same.
  * Note that the @l channels determines the number of output channels in any situation
  * and is therefore must be equal to the maximum number of channels between input and 
  * coefficient matrices.
  *
  * @author Ricard Marxer
  *
  * @sa Filter
  */
class Filter {
protected:
  // Internal parameters
  int _channels; 
  int _length;

  // Internal variables
  MatrixXR _ina;
  MatrixXR _inb;

  MatrixXR _a;
  MatrixXR _b;

  MatrixXR _z;
  MatrixXR _samples;

  void setupState();
  void setupCoeffs();

public:
  /**
     Constructs a band pass filter object with the given @a channels, @a b,
     and @a a coefficients given.
  */
  Filter(int channels = 1);
  Filter(const MatrixXR& b, const MatrixXR& a, int channels);

  void setup();
  void reset();

  /**
     Performs a filtering of each of the columns of @a samples.
     Puts the resulting filtered in the columns of @a filtered.
     
     @param samples matrix of Real values.  A column represents a channel and a row 
     represents a time index.
     
     @param filtered pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows and columns as @a samples.
     
     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
     
     @sa IIRFilter::process
  */
  void process(const MatrixXR& samples, MatrixXR* filtered);

  /**
     @property Filter::a
     @brief the a coefficients of the filter
     
     Both the b and a coefficients matrices are normalized
     by the first row of the a coefficients matrix.
     
     Note that if the first row of the a coefficients matrix
     has elements to zero, some of the filtered samples will 
     result in NaN.
     
     Note that the number of columns in a and b must be the same,
     and that it must be equal to 1 or Filter::channels.

     By default it is a single element matrix of value 1.0.

     @sa Filter::b
  */
  void a( MatrixXR* a ) const;
  void setA( const MatrixXR& a, bool callSetup = true );

  /**
     @property Filter::b
     @brief the b coefficients of the filter
     
     Both the b and a coefficients matrices are normalized
     by the first row of the a coefficients matrix.
     
     Note that if the first row of the a coefficients matrix
     has elements to zero, some of the filtered samples will 
     result in NaN.

     Note that the number of columns in a and b must be the same,
     and that it must be equal to 1 or Filter::channels.
     
     By default it is a single element matrix of value 1.0.

     @sa Filter::a
  */
  void b( MatrixXR* b ) const;
  void setB( const MatrixXR& b, bool callSetup = true );

  /**
     @property Filter::channels
     @brief the number of output channles of the filter
     
     Note that the number of channels must be equal to the 
     number of columns in a and b or to the number of columns
     in the input matrix.
     
     By default it is 1.
  */
  int channels() const;
  void setChannels( int channels, bool callSetup = true );


  /**
     Returns the length of the filter, which is 
     the maximum number of rows between the a and b coefficients 
     matrices.
  */
  int length() const;
};

#endif  /* FILTER_H */
