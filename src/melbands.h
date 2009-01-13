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

#include "typedefs.h"
#include "debug.h"

#include <Eigen/StdVector>

#include "bands.h"

class MelBands {
public:
  Real _lowFreq;
  Real _highFreq;
  int _numBands;
  Real _samplerate;
  int _fftLength;
  
  Bands _bands;

  void triangleWindow(MatrixXR* window, Real start, Real stop, Real center = -1, Real height = Real(1.0));

  Real linearToMelGreenwood1990(Real linearFreq);
  Real melToLinearGreenwood1990(Real melFreq);
  void linearToMelMatrixGreenwood1990(MatrixXR linearFreq, MatrixXR* melFreq);
  void melToLinearMatrixGreenwood1990(MatrixXR melFreq, MatrixXR* linearFreq);

  Real linearToMelStevens1937(Real linearFreq);
  Real melToLinearStevens1937(Real melFreq);
  void linearToMelMatrixStevens1937(MatrixXR linearFreq, MatrixXR* melFreq);
  void melToLinearMatrixStevens1937(MatrixXR melFreq, MatrixXR* linearFreq);

  Real linearToMelFant1968(Real linearFreq);
  Real melToLinearFant1968(Real melFreq);
  void linearToMelMatrixFant1968(MatrixXR linearFreq, MatrixXR* melFreq);
  void melToLinearMatrixFant1968(MatrixXR melFreq, MatrixXR* linearFreq);


public:
  MelBands(Real lowFreq, Real highFreq, int numBands, Real samplerate, int fftLength);

  void setup();

  void process(MatrixXR spectrum, MatrixXR* bands);
  
  void reset();

  std::vector<MatrixXR> weights() const;

  void bandWeights(int band, MatrixXR* bandWeights) const;

  void starts(MatrixXI* result) const;

  int bands() const;
};

#endif  /* MELBANDS_H */
