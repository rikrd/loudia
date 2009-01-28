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

#ifndef LPC_H
#define LPC_H

#include "typedefs.h"
#include "debug.h"

class LPC {
protected:
  // Internal parameters
  int _frameSize;
  int _numCoeffs;
  
  // Internal variables
  MatrixXR _temp;
  MatrixXR _acorr;
  MatrixXR _fullAcorr;

public:
  LPC(int frameSize, int numCoeffs);

  ~LPC();

  void setup();

  void process(const MatrixXR& frame, MatrixXR* lpcCoeffs, MatrixXR* reflectionCoeffs, MatrixXR* error);

  void reset();

  int numCoeffs() const;
};

#endif  /* LPC_H */
