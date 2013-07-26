/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 3 of the License, or   
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

#include "Typedefs.h"
#include "Debug.h"

#include "LPC.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

LPC::LPC(int inputSize, int coefficientCount, Real preEmphasis) 
{
  LOUDIA_DEBUG("LPC: Constructor inputSize: " << inputSize 
        << ", coefficientCount: " << coefficientCount
        << ", preEmphasis: " << preEmphasis);

  if ( coefficientCount > inputSize ) {
    // Thorw ValueError, the number of coefficients must be smaller or equal than the frame size.
  }
  
  setInputSize( inputSize, false ); 
  setCoefficientCount( coefficientCount, false );
  setPreEmphasis( preEmphasis, false );
  
  setup();
}

LPC::~LPC() {}


void LPC::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("LPC: Setting up...");
  
  if ( _preEmphasis != 0.0 ) {
    MatrixXR preCoeffs(2, 1);
    preCoeffs << 1, -_preEmphasis;
    _preFilter.setB( preCoeffs );
  }
  
  _acorrelation.setInputSize( _inputSize, false );
  _acorrelation.setMaxLag( _coefficientCount + 1, false );
  _acorrelation.setUseFft( true, false );
  _acorrelation.setup();
  
  reset();
  
  LOUDIA_DEBUG("LPC: Finished set up...");
}


void LPC::process(const MatrixXR& frame, MatrixXR* lpcCoeffs, MatrixXR* reflectionCoeffs, MatrixXR* error){
  LOUDIA_DEBUG("LPC: Processing...");
  const int rows = frame.rows();
  const int cols = frame.cols();

  if ( cols != _inputSize ) {
    // Throw ValueError, the frames passed are the wrong size
  }
  
  _pre.resize(rows, cols);
  
  if ( _preEmphasis != 0.0 ) {
    for ( int row = 0; row < rows; row++) {
      _preFilter.process( frame.transpose(), &_preRow );
      _pre.row( row ) = _preRow.transpose();
    }
  } else {
    _pre = frame;
  }
  
  LOUDIA_DEBUG("LPC: Processing autocorrelation");
  
  _acorrelation.process(_pre, &_acorr);
  
  LOUDIA_DEBUG("LPC: Processing Levinson-Durbin recursion");

  (*lpcCoeffs).resize(rows, _coefficientCount);
  (*reflectionCoeffs).resize(rows, _coefficientCount - 1);
  (*error).resize(rows, 1);
  
  // Initial values of the LPC coefficients
  (*lpcCoeffs).setZero();
  (*lpcCoeffs).col(0).setOnes();
  
  // Initial value of the Error
  (*error).col(0) = _acorr.col(0);

  (*reflectionCoeffs).setZero();

  for ( int row = 0; row < rows; row++) {  
    Real gamma;
    
    if ((_acorr.array() == 0.).all())
      continue;

    for ( int i = 1; i < _coefficientCount; i++ ) {
      gamma = _acorr(row, i);

      // Use the Eigen reverse()      
//       if ( i >= 2) {
//         gamma += ((*lpcCoeffs).row(row).segment(1, i-1) * _acorr.row(row).segment(1, i-1).transpose().reverse())(0,0);
//       }
      
      // instead of manually walking it in reverse order
      for (int j = 1; j <= i-1; ++j) {
       gamma += (*lpcCoeffs)(row, j) * _acorr(row, i-j);  
      }
      
      // Get the reflection coefficient
      (*reflectionCoeffs)(row, i-1) = - gamma / (*error)(row, 0);

      // Update the error      
      (*error)(row, 0) *= (1 - (*reflectionCoeffs)(row, i-1) * (*reflectionCoeffs).conjugate()(row, i-1));
      
      // Update the LPC coefficients
      if(i >= 2){
        _temp = (*lpcCoeffs).row(row).segment(1, i-1);
        reverseCols(&_temp);
        
        (*lpcCoeffs).row(row).segment(1, i-1) += (*reflectionCoeffs)(row, i-1) * _temp.conjugate();
      }
      
      (*lpcCoeffs)(row, i) = (*reflectionCoeffs)(row, i-1);
    }
  }
  
  LOUDIA_DEBUG("LPC: Finished Processing");
}

void LPC::reset(){
  // Initial values

  if ( _preEmphasis != 0.0 ) {
    _preFilter.reset( );
  }

  _acorrelation.reset( );


}

int LPC::inputSize() const {
  return _inputSize;
}
  
void LPC::setInputSize( int size, bool callSetup ) {
  _inputSize = size;
  if ( callSetup ) setup();
}


int LPC::coefficientCount() const {
  return _coefficientCount;
}

void LPC::setCoefficientCount( int count, bool callSetup ) {
  _coefficientCount = count;
  if ( callSetup ) setup();
}

Real LPC::preEmphasis() const {
  return _preEmphasis;
}

void LPC::setPreEmphasis( Real coefficient, bool callSetup ) {
  _preEmphasis = coefficient;
  if ( callSetup ) setup();
}
