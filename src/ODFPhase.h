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

#ifndef ODFPHASE_H
#define ODFPHASE_H

#include "Typedefs.h"
#include "Debug.h"

#include "ODFBase.h"
#include "Unwrap.h"

class ODFPhase : public ODFBase {
protected:
  // Internal parameters
  int _fftLength;
  bool _weighted;
  bool _normalize;
  
  // Internal variables
  Unwrap _unwrap;

  MatrixXC _spectrum;
  MatrixXR _unwrappedAngle;
  MatrixXR _phaseDiff;
  MatrixXR _instFreq;
  
  void phaseDeviation(const MatrixXC& spectrum, const MatrixXR& spectrumArg, MatrixXR* odfValue);

public:
  ODFPhase(int fftLength, bool weighted = false, bool normalize = false);

  ~ODFPhase();

  void setup();

  void process(const MatrixXC& fft, MatrixXR* odfValue);

  void reset();

};

#endif  /* ODFPHASE_H */
