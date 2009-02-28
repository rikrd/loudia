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

#ifndef SPECTRALNOISESUPPRESSION_H
#define SPECTRALNOISESUPPRESSION_H

#include "Typedefs.h"
#include "Debug.h"

#include "Bands.h"

class SpectralNoiseSuppression {
protected:
  int _fftSize;
  Real _samplerate;
  
  int _k0;
  int _k1;

  MatrixXR _g;
  MatrixXR _noise;

  Bands _bands;

public:
  SpectralNoiseSuppression(int fftSize, Real f0, Real f1, Real samplerate = 1.0);

  ~SpectralNoiseSuppression();

  void setup();

  void process(const MatrixXR& spectrum, MatrixXR* result);

  void reset();
};

#endif  /* SPECTRALNOISESUPPRESSION_H */
