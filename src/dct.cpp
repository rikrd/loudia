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

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>
#include <cmath>
#include "dct.h"

#include "debug.h"
#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

DCT::DCT(int inputLength, int dctLength, bool scale) {
  DEBUG("DCT: Construction inputLength: " << inputLength << ", dctLength: " << dctLength);
  
  if (inputLength < dctLength) {
    // TODO: Throw an exception since dctLength is the number of coefficients to output and it cannot output more
    return;
  }

  _inputLength = inputLength;
  _dctLength = dctLength;
  _scale = scale;
}

DCT::~DCT(){}

void DCT::setup(){
  // Prepare the buffers
  DEBUG("DCT: Setting up...");
  _dctMatrix.resize(_inputLength, _inputLength);

  Real norm = 1.0;
  if ( _scale ) norm = sqrt(Real(2.0)/Real(_inputLength));

  for(int i=0; i < _dctMatrix.rows(); i++){
    for(int j=0; j < _dctMatrix.cols(); j++){
      _dctMatrix(i,j) = norm * cos(Real(i) * M_PI / Real(_inputLength) * (Real(j) + Real(0.5)));
    }
  }
  
  reset();
  DEBUG("DCT: Finished setup.");
}


void DCT::process(MatrixXR input, MatrixXR* dctCoeffs){
  DEBUG("DCT: Processing input.rows(): " << input.rows() << ", dctCoeffs.rows(): " << (*dctCoeffs).rows());
  for ( int i = 0 ; i < input.rows(); i++) {
    // FIX: remove the need of two transpose() calls by calculating the DCT transformation matrix correctly at the beginning
    (*dctCoeffs).row(i) = (_dctMatrix * input.row(i).transpose()).block(0, 0, _dctLength, 1).transpose();
  }
}

void DCT::reset(){
  // Initial values
}
