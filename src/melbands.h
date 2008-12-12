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

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

#include "spectralbands.h"

#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class MelBands: public SpectralBands {
protected:
  Real _lowFreq;
  Real _highFreq;
  int _numBands;
  Real _samplerate;
  int _spectrumLength;

  void triangleWindow(MatrixXR* window, Real start, Real stop, Real center = -1, Real height = Real(1.0));

  Real linearToMelRealStevens1937(Real linearFreq);
  Real melToLinearRealStevens1937(Real melFreq);
  MatrixXR linearToMelStevens1937(MatrixXR linearFreq);
  MatrixXR melToLinearStevens1937(MatrixXR melFreq);

  Real linearToMelRealFant1968(Real linearFreq);
  Real melToLinearRealFant1968(Real melFreq);
  MatrixXR linearToMelFant1968(MatrixXR linearFreq);
  MatrixXR melToLinearFant1968(MatrixXR melFreq);


public:
  MelBands(Real lowFreq, Real highFreq, int numBands, Real samplerate, int spectrumLength);

  void setup();
};

#endif  /* MELBANDS_H */
