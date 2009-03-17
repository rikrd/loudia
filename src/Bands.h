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

#ifndef BANDS_H
#define BANDS_H

#include <Eigen/StdVector>
#include <vector>

#include "Typedefs.h"
#include "Debug.h"

/**
  * @class Bands
  *
  * @brief Algorithm to calculate the sum of values in a given set of weighted bands.
  *
  * This class represents an object to perform calculatinos of band values.
  * The value of a band corresponds to the sum of the values of the input array in the band's
  * positions multiplied with the band's weights.
  *
  * In this implementation the positions of a given band are defined by the index of the first 
  * array cell of the band and the size of the weights array.
  * The full configuration of the bands algorithm is defined using a single column matrix for the 
  * starts of the bands and a vector of single row matrices for the weights each band.
  *
  * Note that the number of rows of the starts matrix and the size of the vector of weights must
  * be the same, and this will be the number of bands.
  *
  * @author Ricard Marxer
  *
  * @sa MelBands
  */
class Bands {
protected:
  // Internal parameters
  MatrixXI _starts;
  std::vector<MatrixXR> _weights;

  // Internal variables

public:
  /**
     Constructs a Bands object with the a single band covering the entire array.
  */
  Bands();
  
  /**
     Constructs a Bands object with the specified @a starts and @a
     weights setting.
     
     @param starts single column matrix of Integers that determine the
     first array cell of each band.
     
     @param weights vector of single row matrices of Reals that determine the
     values of the weights of each band.
  */
  Bands(MatrixXI starts, std::vector<MatrixXR> weights);

  /**
     Destroys the Bands algorithm and frees its resources.
  */
  ~Bands();

  /**
     Calculates the bands of @a frames using the specified starts and weights properties.
     
     @param frames matrix of Real values.
     
     @param bands pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows as @a frames and as many columns as the number of bands
     (rows in the starts matrix and elements in the weights vector). 

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR&  frames, MatrixXR* bands);

  void setup();
  void reset();

  /**
     Return the vector of weights.
  */
  std::vector<MatrixXR> weights() const;

  /**
     Return in @a bandWeights the weights of the band given by the index @a band.
  */
  void bandWeights(int band, MatrixXR* bandWeights) const;

  /**
     Return in @a result the single column matrix of start indices of the bands.
  */
  void starts(MatrixXI* result) const;

  /**
     Return number of bands.
  */
  int bands() const;

  /**
     Determines the @a starts positions and @a weights of the bands.
     Note that the number of rows of starts and the size of weights should be the same and
     will determine the number of bands.
  */
  void setStartsWeights(const MatrixXI& starts, std::vector<MatrixXR> weights);
};

#endif  /* BANDS_H */
