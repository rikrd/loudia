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

#include "typedefs.h"
#include "debug.h"

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>
#include <cmath>
#include "dct.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

DCT::DCT(int inputLength, int dctLength, bool scale, DCTType dctType) {
  DEBUG("DCT: Construction inputLength: " << inputLength << ", dctLength: " << dctLength << ", dctType: " << dctType);
  
  if (inputLength < dctLength) {
    // TODO: Throw an exception since dctLength is the number of coefficients to output and it cannot output more
    return;
  }

  _inputLength = inputLength;
  _dctLength = dctLength;
  _scale = scale;

  _dctType = dctType;

  setup();
}

DCT::~DCT(){}

void DCT::setup(){
  // Prepare the buffers
  DEBUG("DCT: Setting up...");
  _dctMatrix.resize(_inputLength, _inputLength);


  switch(_dctType) {
  case I:
    type1Matrix( &_dctMatrix );
    break;

  case II:
    type2Matrix( &_dctMatrix );
    break;

  case OCTAVE:
    typeOctaveMatrix( &_dctMatrix );
    break;

  }

  
  reset();
  DEBUG("DCT: Finished setup.");
}

void DCT::type1Matrix(MatrixXR* dctMatrix) {
  int length = (*dctMatrix).rows();

  Real norm = 1.0;
  if ( _scale ) norm = sqrt(Real(2.0)/Real(length - 1));
  
  for(int i=0; i < length; i++){
    (*dctMatrix)(i, length - 1) = norm * 0.5 * pow(-1, i);
    for(int j=1; j < length-1; j++){
      (*dctMatrix)(i,j) = norm * cos(Real(j * i) * M_PI / Real(length - 1));
    }
  }

  (*dctMatrix).col(0).setConstant(norm * 0.5);

}

void DCT::type2Matrix(MatrixXR* dctMatrix) {
  int length = (*dctMatrix).rows();

  Real norm = 1.0;
  if ( _scale ) norm = sqrt(Real(2.0)/Real(length));
  
  for(int i=0; i < length; i++){
    for(int j=0; j < length; j++){
      (*dctMatrix)(i,j) = norm * cos(Real(j) * M_PI / Real(length) * (Real(i) + 0.5));
    }
  }
}

void DCT::typeOctaveMatrix(MatrixXR* dctMatrix) {
  int length = (*dctMatrix).rows();

  Real norm = 1.0;
  if ( _scale ) norm = sqrt(2.0/Real(length));
  
  for(int i=0; i < length; i++){
    for(int j=1; j < length; j++){
      (*dctMatrix)(i,j) = norm * cos(Real(j) * M_PI / Real(2 * length) * (Real(2 * i - 1)));
    }
  }
  
  (*dctMatrix).col(0).setConstant( norm * sqrt(0.5) );
}

void DCT::process(MatrixXR input, MatrixXR* dctCoeffs){
  (*dctCoeffs).resize(input.rows(), _dctLength);
  
  for ( int i = 0 ; i < input.rows(); i++) {
    (*dctCoeffs).row(i) = (input.row(i) * _dctMatrix).block(0, 0, 1, _dctLength);
  }
}

void DCT::reset(){
  // Initial values
}
