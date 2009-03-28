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

#ifndef FFT_H
#define FFT_H

#include "Typedefs.h"
#include "Debug.h"

#include <fftw3.h>

/**
  * @class FFT
  *
  * @brief Algorithm to perform a Fast Fourier Transform of a vector of Real values.
  *
  * This class represents an object to perform Fast Fourier Transforms (FFT) on Real data.
  * The FFT is a fast implementation of a Discrete Fourier Transform (DFT).
  * The algorithm takes as input N-point vectors of Real values (N being the frame size)
  * and returns (M / 2 + 1) point vectors of Complex values (M being the FFT size).
  *
  * Note that the algorithm works fastest when M is a power of 2.
  *
  * Since the input is Real valued, the FFT will be symmetric
  * and only half of the output is needed.
  * This processing unit will only return the (M / 2 + 1)-point array 
  * corresponding to positive frequencies of the FFT.
  *
  * When M is different than N the input data is zero padded at the end.
  * Alternatively the algorithm can perform an N/2 rotation and zero pad the center
  * before the FFT to allow a zero phase transform.
  * This is done by using the setZeroPhase() method.
  *
  * @author Ricard Marxer
  *
  * @sa FFTComplex, IFFT, IFFTComplex
  */
class FFT{
protected:
  int _fftSize;
  bool _zeroPhase;

  int _halfSize;
  
  Real* _in;
  fftwf_complex* _out;
  fftwf_plan _fftplan;
  

public:
  /**
     Constructs an FFT object with the specified @a fftSize and @a
     zeroPhase setting.
     
     @param fftSize size of the FFT transform must be > 0, 
     it is the target size of the transform.
     The algorithm performs faster for sizes which are a power of 2.
     
     @param zeroPhase determines whether
     or not to perform the zero phase transform.
  */
  FFT(int fftSize = 1024, bool zeroPhase = true);
  
  /**
     Destroys the FFT algorithm and frees its resources.
  */
  ~FFT();
  
  /**
     Performs a Fast Fourier Transform on each of the rows of @a frames and
     puts the resulting FFT in the rows of @a fft.
     
     @param frames matrix of Real values.  The number of columns of @a frames must
     be smaller or equal to the fftSize property.
     
     @param fft pointer to a matrix of Complex values for the output.  The matrix should
     have the same number of rows as @a frames and (fftSize / 2) + 1 columns. 

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& frames, MatrixXC* fft);

  void setup();
  void reset();

  /**
     Returns the size of the FFT to be performed.  The default is 1024.
     
     @sa setFftSize()
  */
  int fftSize() const;

  /**
     Specifies the @a size of the FFT to be performed.
     The given @a size must be higher than 0.
     Note that if @a size is a power of 2 will perform faster.
     
     @sa fftSize()
  */
  void setFftSize( int size, bool callSetup = true );

  /**
     Returns the zero phase setting.  The default is True.
     
     @sa setZeroPhase()
  */
  bool zeroPhase() const;

  /**
     Specifies the @a zeroPhase setting.
     
     @sa zeroPhase()
  */
  void setZeroPhase( bool zeroPhase, bool callSetup = true );

};

#endif  /* FFT_H */
