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

#ifndef MELBANDS_H
#define MELBANDS_H

#include <Eigen/StdVector>

#include "Typedefs.h"
#include "Debug.h"

#include "Bands.h"

class MelBands {
public:
  enum ScaleType {
    STEVENS = 0,
    FANT = 1,
    GREENWOOD = 2
  };

protected:
  Real _lowFreq;
  Real _highFreq;
  int _numBands;
  Real _samplerate;
  int _fftLength;
  ScaleType _scaleType;

  Bands _bands;

  Real (*_linearToMel)(Real linearFreq);
  
  Real (*_melToLinear)(Real melFreq);
  
  void (*_linearToMelMatrix)(const MatrixXR& linearFreq, MatrixXR* melFreq);
  
  void (*_melToLinearMatrix)(const MatrixXR& melFreq, MatrixXR* linearFreq);

  void triangleWindow(MatrixXR* window, Real start, Real stop, Real center = -1, Real height = Real(1.0));
  
public:
  MelBands(Real lowFreq, Real highFreq, int numBands, Real samplerate, int fftLength, ScaleType scaleType = GREENWOOD);

  void setup();

  void process(const MatrixXR& spectrum, MatrixXR* bands);
  
  void reset();

  std::vector<MatrixXR> weights() const;

  void bandWeights(int band, MatrixXR* bandWeights) const;

  void starts(MatrixXI* result) const;

  int bands() const;

};

#endif  /* MELBANDS_H */
