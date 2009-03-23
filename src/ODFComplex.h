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

#ifndef ODFCOMPLEX_H
#define ODFCOMPLEX_H

#include "Typedefs.h"
#include "Debug.h"

#include "ODFBase.h"
#include "Unwrap.h"

class ODFComplex : public ODFBase {
protected:
  // Internal parameters
  int _fftSize;
  int _halfSize;
  bool _rectified;
  
  // Internal variables
  Unwrap _unwrap;

  MatrixXC _unwrappedSpectrum;
  MatrixXR _unwrappedAngle;
  MatrixXC _spectrumPredict;
  MatrixXR _predictionError;
  
  void spectralDistanceEuclidean(const MatrixXC& spectrum, const MatrixXR& spectrumAbs, const MatrixXR& spectrumArg, MatrixXR* odfValues);

public:
  ODFComplex(int fftSize, bool rectified = false);

  ~ODFComplex();

  void setup();
  void reset();

  void process(const MatrixXC& fft, MatrixXR* odfValue);


};

#endif  /* ODFCOMPLEX_H */
