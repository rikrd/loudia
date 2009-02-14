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

class IFFT{
protected:
  int _fftSize;
  int _frameSize;
  bool _zeroPhase;

  int _halfSize;

  fftwf_complex* _in;
  Real* _out;

  fftwf_plan _fftplan;
  

public:
  IFFT(int fftSize, int frameSize, bool zeroPhase = true);
  ~IFFT();
  
  void process(const MatrixXC& fft, MatrixXR* frame);
  
  void setup();
  void reset();

  int frameSize() const;
  int fftSize() const;
};

#endif  /* IFFT_H */
