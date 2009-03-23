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

#ifndef DCT_H
#define DCT_H

#include "Typedefs.h"
#include "Debug.h"

/**
  * @class DCT
  *
  * @brief Algorithm to perform a Discrete Cosine Transform of a vector of Real values.
  *
  * This class represents an object to perform Discrete Cosine Transform (DCT) on Real data.
  * The algorithm takes as input N-point vectors of Real values 
  * and returns M-point vectors of Real values.
  *
  * 5 types of DCT are implemented:
  * -# Type I
  * -# Type II
  * -# Type III
  * -# Type IV
  * -# Octave's Implementation
  *
  * The DCT type can be selected using the 
  * setDCTType() taking as argument a DCTType.
  *
  *
  * @author Ricard Marxer
  *
  * @sa FFT
  */
class DCT {
public:
  /**
    @enum DCTType
    @brief Specifies the type of the DCT.

    @sa dctType
  */
  enum DCTType {
    I = 0 /**< DCT Type-I */,
    II = 1 /**< DCT Type-II */,
    III = 2 /**< DCT Type-III */,
    IV = 3 /**< DCT Type-IV */,
    OCTAVE = 4 /**< Octave's implementation */
  };

protected:
  // Internal parameters
  int _inputSize;
  int _dctSize;
  Real _scale;

  DCTType _dctType;

  // Internal variables
  MatrixXR _dctMatrix;

  void type1Matrix(MatrixXR* dctMatrix);

  void type2Matrix(MatrixXR* dctMatrix);

  void typeOctaveMatrix(MatrixXR* dctMatrix);

public:
  /**
     Constructs a DCT object with the given @a inputSize, @a dctSize,
     @a scale, @a dctType parameters
     given.
  */
  DCT(int inputSize = 1024, int dctSize = 1024, bool scale = false, DCTType dctType = OCTAVE);
  
  /**
     Destroys the DCT algorithm and frees its resources.
  */
  ~DCT();

  void reset();
  void setup();

  /**
     Performs a Discrete Cosine Transform on each of the rows of @a frames and
     puts the resulting DCT in the rows of @a dct.
     
     @param frames matrix of Real values.  The number of columns of @a frames must
     be equal to the inputSize property.
     
     @param dct pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows as @a frames and dctSize columns. 

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& frames, MatrixXR* dct);

  /**
     Returns the type of the DCT
     
     By default it is OCTAVE.
  */
  DCTType dctType() const;

  /**
     Specifies the type of the DCT
  */
  void setDctType( DCTType type, bool callSetup = true );

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
     Returns the output size of the DCT
     
     Note that the result will when performing
     the DCT at the most inputSize coefficients
     will be outputed.
     
     By default it is 1024.
  */
  int dctSize() const;

  /**
     Specifies the output size of the DCT
     
     Note that the result will when performing
     the DCT at the most inputSize coefficients
     will be outputed.
     
     By default it is 1024.
  */
  void setDctSize( int size, bool callSetup = true );
};

#endif  /* DCT_H */
