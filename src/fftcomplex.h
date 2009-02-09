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

#ifndef FFTCOMPLEX_H
#define FFTCOMPLEX_H

#include "typedefs.h"
#include "debug.h"

#include <fftw3.h>

class FFTComplex{
protected:
  int _frameSize;
  int _fftSize;
  bool _zeroPhase;

  fftwf_complex* _in;
  
  fftwf_complex* _out;

  fftwf_plan _fftplan;
  
  template <typename FrameMatrixType>
  void process(const FrameMatrixType& frames, MatrixXC* fft);


public:
  FFTComplex(int frameSize, int fftSize, bool zeroPhase = true);
  ~FFTComplex();
  
  void process(const MatrixXC& frames, MatrixXC* fft);
  void process(const MatrixXR& frames, MatrixXC* fft);
  
  void setup();
  void reset();

  int frameSize() const;
  int fftSize() const;
};

#endif  /* FFTCOMPLEX_H */
