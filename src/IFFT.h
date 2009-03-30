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

#ifndef IFFT_H
#define IFFT_H

#include "Typedefs.h"
#include "Debug.h"

#include <fftw3.h>

/**
  * @class IFFT
  *
  * @brief Algorithm to perform an Inverse Fast Fourier Transform of a vector of Complex 
  * values representing the positive frequencies half of a symmetric FFT.
  *
  * The IFFT is a fast implementation of an Inverse Discrete Fourier Transform (IDFT).
  * The algorithm takes as input (M / 2 + 1) point vectors of Complex 
  * values (M being the FFT size), and returns N point vectors of Real 
  * values (N being the frame size).
  *
  * The input of the IFFT is assumed to be the positive frequencies half of an M point
  * magnitude symmetric and phase antisymmetric FFT.  Therefore the result is a Real value 
  * vector.
  *
  * Note that N can be smaller than M.
  * In this case the last ( M - N ) coefficients
  * will be discarded, since it assumes that zero padding has been made
  * at the end of the frame prior to the forward FFT transform.
  *
  * Alternatively the algorithm can undo the center zeropadding and
  * the N/2 rotation if done durnig the FFT forward transform.
  * This is specified by using the setZeroPhase() method.
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
     
     @param zeroPhase specifies whether
     or not the zero phase method was performed.
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
     where fftSize is parameter of the constructor or specified by setFftSize().
     
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
