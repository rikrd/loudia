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

#ifndef ODFCOG_H
#define ODFCOG_H

#include "typedefs.h"
#include "debug.h"

#include "odfbase.h"
#include "peakdetect.h"

class ODFCOG : public ODFBase {
protected:
  // Internal parameters
  int _fftLength;
  int _peakCount;
  int _bandwidth;
  
  // Internal variables
  MatrixXR _spectrumAbs2;
  MatrixXR _spectrumArg;
  MatrixXR _spectrumArgDeriv;

  MatrixXR _peakPos;
  MatrixXR _peakMag;
  MatrixXR _peakArg;

  MatrixXR _cog;

  PeakDetect _peaker;
  
public:
  ODFCOG(int fftLength, int bandwidth = 6, int peakCount = 40);

  ~ODFCOG();

  void setup();

  void process(const MatrixXC& fft, MatrixXR* odfValue);

  void reset();

};

#endif  /* ODFCOG_H */
