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

#ifndef PITCHINVERSEPROBLEM_H
#define PITCHINVERSEPROBLEM_H

#include "Typedefs.h"
#include "Debug.h"

class PitchInverseProblem {
protected:
  int _fftSize;
  int _halfSize;
  Real _f0;
  Real _f1;
  Real _fPrec;
  int _numHarmonics;
  int _numMaxPitches;
  int _peakBandwidth;

  Real _samplerate;

  Real _tMin;
  Real _tMax;
  Real _tPrec;
  Real _alpha;
  Real _beta;
  Real _inharmonicity;

  MatrixXR _projectionMatrix;
  MatrixXR _inverseProjectionMatrix;

public:
  PitchInverseProblem(int fftSize, Real f0, Real f1, Real samplerate = 1.0, Real fPrec = 0.01, int numHarmonics = 5, int numMaxPitches = 5, int peakBandwidth = 8);

  ~PitchInverseProblem();

  void setup();

  void process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies);

  void projectionMatrix(MatrixXR* matrix) const;

  void inverseProjectionMatrix(MatrixXR* inverseMatrix) const;

  Real harmonicWeight(Real period, Real tLow, Real tUp, int harmonicIndex);

  Real harmonicSpread(Real perdio, Real tLow, Real tUp, int harmonicIndex);

  Real harmonicPosition(Real period, Real tLow, Real tUp, int harmonicIndex);

  void reset();
};

#endif  /* PITCHINVERSEPROBLEM_H */
