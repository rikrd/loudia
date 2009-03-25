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

#ifndef SPECTRALODF_H
#define SPECTRALODF_H

#include "Typedefs.h"
#include "Debug.h"

#include "SpectralODFBase.h"

/**
  * @class SpectralODF
  *
  * @brief Algorithm to perform onset detection functions on vectors of Complex values representing the FFT of a signal.
  *
  * This class wraps several kind of spectral Onset Detection Functions (ODF).  A spectral onset detection function is 
  * a function mapping from an FFT of an audio signal to a Real value.  The onset detection function represents a value
  * proportional to the probability of the frame being an onset.
  *
  * The algorithm takes as input N-point vectors of Complex values 
  * and returns Real values.
  *
  * 5 types of ODF methods are implemented:
  * -# Spectral flux
  * -# High frequency content
  * -# Phase deviation
  * -# Weighted phase deviation
  * -# Normalized weighted phase deviation
  * -# Modified Kullback-Liebler
  * -# Complex domain
  * -# Rectified complex domain
  * -# Peak center of gravity
  *
  * The Phase deviation, Weighted phase deviation, Normalized weighted phase deviation, Complex domain and 
  * Rectified complex domain methods require the 2 past FFT frames.
  * Therefore if any of these methods are specified, the process() method would require an input matrix of at least 3 rows 
  * and will output ODF values for all the rows but the first two.
  *
  * The Spectral flux method requires the past FFT frame.
  * Therefore if this method is used, the process() method would require an input matrix of at least 2 rows 
  * and will output ODF values for all the rows but the first.
  * 
  * The ODF method can be selected using the 
  * setOdfMethod() taking as argument an ODFMethod.
  *
  * This function is often use to perform beat estimation and event segmentation tasks.
  *
  * @author Ricard Marxer
  *
  * @sa FFT
  */
class SpectralODF : SpectralODFBase {
public:
  /**
    @enum ODFMethod
    @brief Specifies the method for calculating the onset detection function to be used.

    @sa odfMethod
  */
  enum ODFMethod {
    FLUX = 0                          /**< Spectral flux method  */,
    HIGH_FREQUENCY_CONTENT = 1        /**< High frequency content method  */,
    PHASE_DEVIATION = 2               /**< Phase deviation method  */,
    WEIGHTED_PHASE_DEVIATION = 3      /**< Weighted phase deviation method  */,
    NORM_WEIGHTED_PHASE_DEVIATION = 4 /**< Normalized weighted phase deviation method  */,
    MODIFIED_KULLBACK_LIEBLER = 5     /**< Modified Kullback-Liebler method  */,
    COMPLEX_DOMAIN = 6                /**< Complex domain method  */,
    RECTIFIED_COMPLEX_DOMAIN = 7      /**< Rectified complex domain method  */,
    CENTER_OF_GRAVITY = 8             /**< Peak center of gravity method  */
  };

protected:
  // Internal parameters
  ODFMethod _odfMethod;
  
  // Internal variables
  SpectralODFBase* _odf;

public:
  /**
     Constructs a spectral onset detection function object with the specified @a fftSize and 
     @a odfMethod settings.
     
     @param fftSize size of the input FFT frames, must be > 0.
     
     @param odfMethod the onset detection method to be used
  */
  SpectralODF(int fftSize = 1024, ODFMethod odfMethod = COMPLEX_DOMAIN);
  
  /**
     Destroys the algorithm and frees its resources.
  */
  ~SpectralODF();

  void setup();
  void reset();

  /**
     Calculates the onset detection method on each of the rows of @a ffts and
     puts the resulting onset detection function values in the rows of @a odfs.
     
     @param ffts matrix of Complex values.  The number of columns of @a ffts 
     must be equal to the fftSize / 2 + 1.  Some onset detection methods require
     a minimum of 2 (or 3) rows to calculate the one onset detection function value.
     In this cases onset detection values for the first 2 (or 3 respectively) rows
     will not be output.
     
     @param odfs pointer to a single-column matrix of Real values for the output.  The matrix should
     have the same number of rows as @a ffts (minus 1 or 2 depending on the method used) and 1 single column. 

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXC& ffts, MatrixXR* odfs);
  
  /**
     Returns the method for calculating the spectral onset detection function.
     
     By default it is COMPLEX_DOMAIN.
  */
  ODFMethod odfMethod() const;

  /**
     Specifies the @a method for calculating the spectral onset detection function.
     
     Note that PHASE_DEVIATION, WEIGHTED_PHASE_DEVIATION, NORM_WEIGHTED_PHASE_DEVIATION, 
     COMPLEX_DOMAIN and RECTIFIED_COMPLEX_DOMAIN methods require at least 3 FFT frames, and therefore
     the input matrix to the process() method must have at least 3 rows.  In these cases
     the output matrix will be 2 rows smaller than the input, since it only calculates ODF values
     for all the rows but the first two.

     Note that SPECTRAL_FLUX method requires at least 2 FFT frames, and therefore
     the input matrix to the process() method must have at least 2 rows.  In this case
     the output matrix will be 1 row smaller than the input, since it only calculates ODF values
     for all the rows but the first.
     
     @param method the method used for calculating the spectral onset detection function.
     
     @param callSetup a flag specifying whether the setup() method must be called after setting the parameter.
  */
  void setOdfMethod( ODFMethod method, bool callSetup = true );

};

#endif  /* SPECTRALODF_H */
