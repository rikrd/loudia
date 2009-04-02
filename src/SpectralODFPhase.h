/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 3 of the License, or   
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

#ifndef SPECTRALODFPHASE_H
#define SPECTRALODFPHASE_H

#include "Typedefs.h"
#include "Debug.h"

#include "SpectralODFBase.h"
#include "Unwrap.h"

class SpectralODFPhase : public SpectralODFBase {
protected:
  // Internal parameters
  int _fftSize;
  int _halfSize;
  bool _weighted;
  bool _normalize;
  
  // Internal variables
  Unwrap _unwrap;
  
  MatrixXR _unwrappedAngle;
  MatrixXR _phaseDiff;
  MatrixXR _instFreq;
  
  void phaseDeviation(const MatrixXC& spectrum, const MatrixXR& spectrumArg, MatrixXR* odfValue);

public:
  SpectralODFPhase(int fftSize, bool weighted = false, bool normalize = false);

  ~SpectralODFPhase();

  void setup();

  void process(const MatrixXC& fft, MatrixXR* odfValue);

  void reset();

};

#endif  /* SPECTRALODFPHASE_H */
