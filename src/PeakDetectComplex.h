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

#ifndef PEAKDETECTCOMPLEX_H
#define PEAKDETECTCOMPLEX_H

#include "Typedefs.h"
#include "Debug.h"

#include <limits>

class Peakdetectcomplex {
public:
  enum SortType {
    NOSORT              = 0,
    BYMAGNITUDE         = 1,
    BYPOSITION          = 2
  };

protected:
  // Internal parameters
  int _numPeaks;
  int _minPeakWidth;
  int _numCandidates;
  Real _minPeakContrast;
    
  SortType _sort;

  // Internal variables
  MatrixXR _magnitudes;
  MatrixXR _phases;

public:
  PeakDetectComplex(int numPeaks, SortType sort = BYMAGNITUDE, int minPeakWidth = 3, int numCandidates = -1, Real minPeakContrast = 0);

  ~PeakDetectComplex();

  void setup();

  void process(const MatrixXC& fft,
               MatrixXR* peakPositions, MatrixXR* peakMagnitudes, MatrixXR* peakPhases);

  void reset();

  int numPeaks() const;

  int minPeakWidth() const;

};

#endif  /* PEAKDETECTCOMPLEX_H */
