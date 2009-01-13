/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com.com>
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

#ifndef ODFCOMPLEX_H
#define ODFCOMPLEX_H

#include "typedefs.h"
#include "debug.h"

#include "unwrap.h"

class ODFComplex {
protected:
  // Internal parameters
  int _fftLength;
  
  // Internal variables
  Unwrap _unwrap;

  MatrixXC _spectrum;
  MatrixXC _unwrappedSpectrum;
  MatrixXR _unwrappedAngle;
  MatrixXC _spectrumPredict;
  MatrixXR _predictionError;
  
  Real spectralDistanceEuclidean(MatrixXC spectrum, MatrixXR spectrumAbs, MatrixXR spectrumArg);
  Real spectralDistanceEuclideanWeighted(MatrixXC spectrum, MatrixXR spectrumAbs, MatrixXR spectrumArg);
  Real spectralDistanceHypot(MatrixXC spectrum, MatrixXR spectrumAbs, MatrixXR spectrumArg);

public:
  ODFComplex(int fftLength);

  ~ODFComplex();

  void setup();

  void process(MatrixXC fft, MatrixXR* odfValue);

  void reset();

};

#endif  /* ODFCOMPLEX_H */
