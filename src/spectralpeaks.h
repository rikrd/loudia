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

#ifndef SPECTRALPEAKS_H
#define SPECTRALPEAKS_H

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

#include "melbands.h"
#include "dct.h"

#include "typedefs.h"

//using namespace std;

// import most common Eigen types 
//using namespace Eigen;

class SpectralPeaks {
protected:
  // Internal parameters
  int _numPeaks;
    
  // Internal variables

public:
  SpectralPeaks(int numPeaks);

  ~SpectralPeaks();

  void setup();

  void process(MatrixXR spectrum, MatrixXR* peakMagnitudes, MatrixXi* peakPositions);

  void reset();

  int numPeaks() const;

};

#endif  /* SPECTRALPEAKS_H */
