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

#ifndef PITCHACF_H
#define PITCHACF_H

#include "Typedefs.h"
#include "Debug.h"

#include "PeakDetect.h"
#include "PeakInterpolate.h"
#include "Autocorrelation.h"

class PitchACF {
protected:
  int _fftSize;
  int _halfSize;

  Real _samplerate;

  PeakDetect _peak;
  PeakInterpolate _peakInterp;
  Autocorrelation _acorr;

  MatrixXR _acorred;
  MatrixXC _acorredC;

public:
  PitchACF(int fftSize, Real samplerate = 1.0);

  ~PitchACF();

  void setup();

  void process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies);

  void reset();
};

#endif  /* PITCHACF_H */
