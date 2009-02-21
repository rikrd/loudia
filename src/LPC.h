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

#ifndef LPC_H
#define LPC_H

#include "Typedefs.h"
#include "Debug.h"

#include "Filter.h"

class LPC {
protected:
  // Internal parameters
  int _frameSize;
  int _numCoeffs;
  Real _preEmphasis;
  
  // Internal variables
  MatrixXR _pre;
  MatrixXR _preRow;
  MatrixXR _temp;
  MatrixXR _acorr;
  Filter _preFilter;

public:
  /**
   *
   * The LPC often uses a preemphasis filter to emphasize the
   * higher freqs.  The filter is of the form a = [1 a1]
   *
   * with usually 0.96 <= a1 <= 0.99
   *
   */
  LPC(int frameSize, int numCoeffs, Real preEmphasis = 0.0);

  ~LPC();

  void setup();

  void process(const MatrixXR& frame, MatrixXR* lpcCoeffs, MatrixXR* reflectionCoeffs, MatrixXR* error);

  void reset();

  int numCoeffs() const;
};

#endif  /* LPC_H */
