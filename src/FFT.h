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
  * @brief Fast Fourier Transform processor unit for Real data.
  *
  * This class represents an object to perform Fast Fourier Transforms on Real data.
  * It allows to calculate the Discrete Fourier Transform (DFT) of N-point vectors of Real values,
  * returning M-point vectors of Complex values.
  *
  *
  * The algorithm works fastest when M is a power of 2.
  * When M is different than N the input data is zero padded.
  *
  * @link #zerophase zero phase method@endlink
  * Optionally the algorithm can perform a N/2 rotation of the
  * input data and zero pad it in the center of the rotated data to
  * allow a zero phase DFT transformation.
  *
  * Since the input is Real valued, the DFT will be symmetric
  * and only half of the output is needed.
  * This processing unit is the (M / 2 + 1)-point array corresponding
  * to positive frequencies.
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

    @a fftSize must be > 0, it is the target size of the transform.
    The algorithm performs faster for sizes which are a power of 2.

    The @a zeroPhase setting determines whether
    or not to perform the @l{FFT}{zero phase} method.
  */
  FFT(int fftSize, bool zeroPhase = true);
  ~FFT();
  
  void process(const MatrixXR& frames, MatrixXC* fft);
  
  void setup();
  void reset();
  
  int fftSize() const;
};

#endif  /* FFT_H */
