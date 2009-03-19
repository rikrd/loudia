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

#ifndef IFFT_H
#define IFFT_H

#include "Typedefs.h"
#include "Debug.h"

#include <fftw3.h>

/**
  * @class IFFT
  *
  * @brief Inverse Fast Fourier Transform (IFFT) algorithm for an FFT of Real data.
  *
  * This class represents an object to perform 
  * Inverse Fast Fourier Transforms (IFFT) on an FFT of Real data.
  *
  * The IFFT is a fast implementation of an Inverse Discrete Fourier Transform (IDFT).
  * The algorithm takes as input (M / 2 + 1)-point vectors of Complex values 
  * and returns M-point vectors of Real values.
  *
  * @link #zerophase zero phase method@endlink
  * Optionally the algorithm can perform a N/2 rotation of the
  * input data and zero pad it in the center of the rotated data to
  * allow a zero phase DFT transformation.
  *
  * @author Ricard Marxer
  *
  * @sa FFT, FFTComplex, IFFTComplex
  */
class IFFT{
protected:
  int _fftSize;
  bool _zeroPhase;

  int _halfSize;

  fftwf_complex* _in;
  Real* _out;

  fftwf_plan _fftplan;
  

public:
  /**
     Constructs an IFFT object with the specified @a fftSize and @a
     zeroPhase setting.
     
     @param fftSize size of the IFFT transform must be > 0, 
     it is the target size of the transform.
     The algorithm performs faster for sizes which are a power of 2.
     
     @param zeroPhase determines whether
     or not the @l{FFT}{zerophase} method was performed.
  */
  IFFT(int fftSize = 1024, bool zeroPhase = true);
  
  /**
     Destroys the IFFT algorithm and frees its resources.
  */
  ~IFFT();
  
  /**
     Performs a Inverse Fast Fourier Transform on each of the rows of @a fft and
     puts the resulting IFFT in the rows of @a frames.
     
     @param fft matrix of Complex values.  The number of columns of @a fft must
     be equal to the (fftSize / 2) + 1, 
     where fftSize is parameter of the constructor or specified by IFFT::setFftSize().
     
     @param frame pointer to a matrix of Real values for the output.  The matrix should
     have the same number of rows as @a fft and fftSize columns.  

     Note that if the zeroPhase setting is true, the resulting IFFT transforms
     will be rotated to compensate for Zero Phase method that may have been performed
     when the FFT had been done.

     Note that if the output matrix is not of the required size it will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXC& fft, MatrixXR* frame);
  
  void setup();
  void reset();

  /**
     Returns the size of the FFT to be performed.  The default is 1024.
     
     @sa setFftSize()
  */
  int fftSize() const;

  /**
     Specifies the @a size of the IFFT to be performed.
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

#endif  /* IFFT_H */
