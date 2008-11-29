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

#ifndef DCT_H
#define DCT_H

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>

#include "typedefs.h"

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class DCT {
protected:
  // Internal parameters
  int _inputLength;
  int _dctLength;
  Real _scale;

  // Internal variables
  MatrixXR _dctMatrix;
  
public:
  DCT(int inputLength, int dctLength, bool scale = true);

  ~DCT();

  void setup();

  void process(MatrixXR input, MatrixXR* dctCoeffs);

  void reset();
};

#endif  /* DCT_H */
