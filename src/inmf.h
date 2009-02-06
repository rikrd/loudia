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

#ifndef INMF_H
#define INMF_H

#include "typedefs.h"
#include "debug.h"

class INMF {
protected:
  // Internal parameters
  int _fftSize;
  int _numComponents;

  int _maxIterations;
  Real _maxError;

  Real _eps;

  int _numPast;
  int _numNew;

  // Coefficients between 0 and 1 which represent
  // how much should the past and the new be taken
  // into account
  Real _pastCoeff;
  Real _newCoeff;
  
  // Internal variables
  MatrixXR _H, _V, _W, _VH, _HH;

public:
  INMF(int fftSize, int numComponents, int numPast, Real pastCoeff, Real newCoeff,  int maxIterations = 10, Real maxError = 10, Real eps = 1e-9);

  ~INMF();

  void setup();

  void process(const MatrixXR& v, MatrixXR* w, MatrixXR* h);

  void reset();

};

#endif  /* NMF_H */
