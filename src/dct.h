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

#include "typedefs.h"
#include "debug.h"

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>



// import most common Eigen types 
//using namespace Eigen;

class DCT {
public:
  enum DCTType {
    I = 0,
    II = 1,
    III = 2,
    IV = 3,
    OCTAVE = 4
  };

protected:
  // Internal parameters
  int _inputLength;
  int _dctLength;
  Real _scale;

  DCTType _dctType;

  // Internal variables
  MatrixXR _dctMatrix;

  void type1Matrix(MatrixXR* dctMatrix);

  void type2Matrix(MatrixXR* dctMatrix);

  void typeOctaveMatrix(MatrixXR* dctMatrix);

public:
  DCT(int inputLength, int dctLength, bool scale = false, DCTType dctType = II);

  ~DCT();

  void setup();

  void process(MatrixXR input, MatrixXR* dctCoeffs);

  void reset();
};

#endif  /* DCT_H */
