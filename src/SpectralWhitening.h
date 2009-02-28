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

#ifndef SPECTRALWHITENING_H
#define SPECTRALWHITENING_H

#include "Typedefs.h"
#include "Debug.h"

#include "MelBands.h"
#include "Resample.h"

class SpectralWhitening {
protected:
  int _fftSize;
  int _halfSize;
  Real _f0;
  Real _f1;

  Real _samplerate;
  Real _compressionFactor;
  int _numBands;

  MelBands::ScaleType _scaleType;

  MatrixXR _centers;

  MatrixXR _bandEnergy;
  MatrixXR _compressionWeights;

  MelBands _bands;

public:
  SpectralWhitening(int fftSize, Real f0, Real f1, Real samplerate = 1.0, Real compressionFactor = 0.33, int numBands = 30, MelBands::ScaleType scaleType = MelBands::GREENWOOD);

  ~SpectralWhitening();

  void setup();

  void process(const MatrixXR& spectrum, MatrixXR* result);

  void reset();
};

#endif  /* SPECTRALWHITENING_H */
