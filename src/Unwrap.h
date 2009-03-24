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

#ifndef UNWRAP_H
#define UNWRAP_H

#include "Typedefs.h"
#include "Debug.h"

/**
  * @class Unwrap
  *
  * @brief Algorithm to unwrap phases vectors represented as vectors of Real values.
  *
  * This class represents an object to unwrap vectors of phases.
  * The algorithm takes as input N-point vectors of Real values 
  * and returns N-point vectors of Real values.
  *
  * Unwrapping consists in removing phase jumps larger than Pi or smaller to -Pi.
  *
  * The size of the input vectors can be specified using setInputSize().
  *
  * @author Ricard Marxer
  *
  * @sa FFT
  */
class Unwrap {
protected:
  MatrixXR _diff;
  MatrixXR _upsteps;
  MatrixXR _downsteps;
  MatrixXR _shift;

  int _inputSize;

public:
  /**
     Constructs an unwrap object with the given @a inputSize.
  */
  Unwrap(int inputSize = 1024);

  /**
     Destroys the algorithm and frees its resources.
  */
  ~Unwrap();

  void setup();
  void reset();

  /**
     Performs the unwrapping on each of the rows of @a phases and
     puts the resulting unwrapped phases in the rows of @a unwrapped.
     
     @param phases matrix of Real values.  The number of columns of @a phases must
     be equal to the inputSize.
     
     @param unwrapped pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows as @a phases and inputSize columns. 

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& phases, MatrixXR* unwrapped);

  /**
     Returns the input size of the algorithm.
     
     By default it is 1024.
  */
  int inputSize() const;
  
  /**
     Specifies the input size of the algorithm.
  */
  void setInputSize( int size, bool callSetup = true );

};

#endif  /* UNWRAP_H */
