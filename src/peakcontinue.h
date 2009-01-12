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

#ifndef PEAKCONTINUE_H
#define PEAKCONTINUE_H

#include "typedefs.h"
#include "debug.h"

class PeakContinue {
protected:
  // Internal parameters
  int _numTrajectories;
  Real _maxFreqBinChange;

  // Internal variables
  MatrixXR _trajPositions, _trajMagnitudes;
  MatrixXR _pastTrajPositions, _pastTrajMagnitudes;
  
  bool createTrajectory(Real peakPos, Real peakMag,
                        MatrixXR* pastTrajPositions, MatrixXR* pastTrajMagnitudes,
                        MatrixXR* trajPositions, MatrixXR* trajMagnitudes,
                        int row);
    

public:
  PeakContinue(int numTrajectories, Real maxFreqBinChange);

  ~PeakContinue();

  void setup();

  void process(MatrixXC fft, 
               MatrixXR peakPositions, MatrixXR peakMagnitudes,
               MatrixXR* trajPositions, MatrixXR* trajMagnitudes);

  void reset();

};

#endif  /* PEAKCONTINUE_H */
