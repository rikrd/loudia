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

#ifndef IFFTCOMPLEX_H
#define IFFTCOMPLEX_H

#include "Typedefs.h"
#include "Debug.h"

#include <fftw3.h>


/**
  * @class IFFTComplex
  *
  * @brief Algorithm to perform an Inverse Fast Fourier Transform of a vector of Complex 
  * values representing the full FFT.
  *
  * The IFFT is a fast implementation of an Inverse Discrete Fourier Transform (IDFT).
  * The algorithm takes as input M point vectors of Complex 
  * values (M being the FFT size), and returns N point vectors of Real 
  * values (N being the frame size).
  *
  * Note that N can be smaller than M.
  * In this case the last ( M - N ) coefficients
  * will be discarded, since it assumes that zero padding has been made
  * at the end of the frame prior to the forward FFT transfor.
  *
  * Alternatively the algorithm can undo the center zeropadding and
  * the N/2 rotation if done durnig the FFT forward transform.
  * This is specified by using the setZeroPhase() method.
  *
  * @author Ricard Marxer
  *
  * @sa FFT, FFTComplex, IFFT
  */
class IFFTComplex{
protected:
  int _fftSize;
  int _frameSize;
  bool _zeroPhase;

  fftwf_complex* _in;
  fftwf_complex* _out;

  fftwf_plan _fftplan;
  
  template <typename FrameMatrixType>
  void process(const FrameMatrixType& ffts, MatrixXC* frames);

public:
  IFFTComplex(int fftSize, int frameSize, bool zeroPhase = true);
  ~IFFTComplex();
  
  void process(const MatrixXC& ffts, MatrixXC* frames);
  void process(const MatrixXR& ffts, MatrixXC* frames);
  
  void setup();
  void reset();

  /**
     Returns the size of the FFT to be processed.
     The default is 1024.
     
     @sa setFftSize()
  */
  int fftSize() const;

  /**
     Specifies the @a size of the FFT to be processed.
     The given @a size must be higher than 0.
     Note that if @a size is a power of 2 will perform faster.
     
     @sa fftSize()
  */
  void setFftSize( int size, bool callSetup = true );

  /**
     Returns the size of the target frame.
     The default is 1024.
     
     @sa setFrameSize()
  */
  int frameSize() const;

  /**
     Specifies the @a size of the target frame.
     The given @a size must be higher than 0.
     Note that if @a size is a power of 2 will perform faster.
     
     @sa frameSize()
  */
  void setFrameSize( int size, bool callSetup = true );

  /**
     Returns the zero phase setting.
     The default is true.
     
     @sa setZeroPhase()
  */
  bool zeroPhase() const;

  /**
     Specifies the @a zeroPhase setting.
     
     @sa zeroPhase()
  */
  void setZeroPhase( bool zeroPhase, bool callSetup = true );
};

#endif  /* IFFTCOMPLEX_H */
