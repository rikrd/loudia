/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
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

#ifndef RESAMPLE_H
#define RESAMPLE_H

#include "Typedefs.h"
#include "Debug.h"

#include <samplerate.h>

/**
  * @class Resample
  *
  * @brief Algorithm to resample vectors represented as vectors of Real values.
  *
  * This class represents an object to resample vectors of Real values.
  * The algorithm takes as input N-point vectors of Real values 
  * and returns M-point vectors of Real values.
  *
  * This algorithm changes the sampling rate of input vectors, modifying their number
  * of samples.
  *
  * The number of samples in the input and output vectors can be specified using setInputSize()
  * and setOutputSize() respectively.
  *
  * The ratio between the output and input sampling rates can be specified
  * using setResamplingRatio().  Usually the sampling ratio will be chosen
  * close to outputSize / inputSize.
  *
  * @author Ricard Marxer
  *
  * @sa FFT
  */
class Resample{
public:
  /**
    @enum ResamplingMethod
    @brief Specifies the resampling method to be used.

    @sa resamplingMethod
  */
  enum ResamplingMethod {    
    SINC_BEST_QUALITY       = 0 /**< Best quality cardinal sine method  */,
    SINC_MEDIUM_QUALITY     = 1 /**< Medium quality cardinal sine method */,
    SINC_FASTEST            = 2 /**< Fastest cardinal sine method */,
    ZERO_ORDER_HOLD         = 3 /**< Hold the value until the next sample */,
    LINEAR                  = 4 /**< Linear interpolation between samples */
  };
  
protected:
  int _inputSize;
  int _outputSize;
  Real _resamplingRatio;
  ResamplingMethod _resamplingMethod;

  SRC_DATA _resampleData;
  
public:
  /**
     Constructs a Resampler object with the specified @a inputSize, 
     @a outputSize, @a resamplingRatio and @a resamplingMethod settings.
     
     @param inputSize size of the input frames to be resampled,
     must be > 0.

     @param outputSize size of the output resampled frames,
     must be > 0.

     @param resamplingRatio the ratio between the output sampling rate 
     and the input sampling rate
     
     @param resamplingMethod the resampling method to be used
  */
  Resample(int inputSize = 1024, int outputSize = 1024, Real resamplingRatio = 1.0, ResamplingMethod resamplingMethod = SINC_BEST_QUALITY);

  /**
     Destroys the algorithm and frees its resources.
  */  
  ~Resample();

  /**
     Performs the resampling of each of the rows of @a frames and
     puts the resulting resampled frames in the rows of @a resampled.
     
     @param frames matrix of Real values.  The number of columns of @a frames 
     must be equal to the inputSize.
     
     @param resampled pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows as @a frames and outputSize columns. 

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */  
  void process(const MatrixXR& frames, MatrixXR* resampled);
  
  void setup();
  void reset();

  /**
     Returns the input size of the algorithm.
     
     By default it is 1024.
  */
  int inputSize() const;
  
  /**
     Specifies the input @a size of the algorithm.
  */
  void setInputSize( int size, bool callSetup = true );

  /**
     Returns the output size of the algorithm.
     
     By default it is 1024.
  */
  int outputSize() const;
  
  /**
     Specifies the output @a size of the algorithm.
  */
  void setOutputSize( int size, bool callSetup = true );

  /**
     Returns the ratio between the output and input sampling rate.
     Note that this value is normally around outputSize / inputSize.

     By default it is 1.0.
  */
  Real resamplingRatio() const;

  /**
     Specifies the @a ratio between the output and input sampling rate.
  */
  void setResamplingRatio( Real ratio, bool callSetup = true );

  /**
     Returns the resampling method to be used.
     
     By default it is SINC_BEST_QUALITY.
  */
  ResamplingMethod resamplingMethod() const;

  /**
     Specifies the resampling @a method to be used.
  */
  void setResamplingMethod( ResamplingMethod method, bool callSetup = true );

};

#endif  /* RESAMPLE_H */
