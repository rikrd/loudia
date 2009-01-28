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

#include "typedefs.h"
#include "debug.h"

#include "lpc.h"
#include "utils.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

LPC::LPC(int frameSize, int numCoeffs) : 
  _frameSize(frameSize), 
  _numCoeffs(numCoeffs)
{
  DEBUG("LPC: Constructor frameSize: " << frameSize <<
        ", numCoeffs: " << numCoeffs);
  
  setup();
}

LPC::~LPC() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void LPC::setup(){
  // Prepare the buffers
  DEBUG("LPC: Setting up...");

  reset();
  
  DEBUG("LPC: Finished set up...");
}


void LPC::process(const MatrixXR& frame, MatrixXR* lpcCoeffs, MatrixXR* reflectionCoeffs, MatrixXR* error){
  DEBUG("LPC: Processing...");
  const int rows = frame.rows();
  const int cols = frame.cols();
  
  if ( cols != _frameSize ) {
    // Throw ValueError, the frames passed are the wrong size
  }

  DEBUG("LPC: Processing autocorrelation");
  
  correlate(frame, frame, &_fullAcorr);
  
  _acorr.resize(rows, _numCoeffs + 1);
  _acorr.setZero();
  _acorr.block(0, 0, rows, _numCoeffs) = _fullAcorr.block(0, _frameSize - 1, rows, _numCoeffs);

  DEBUG("LPC: Processing Levinson-Durbin recursion");

  (*lpcCoeffs).resize(rows, _numCoeffs);
  (*reflectionCoeffs).resize(rows, _numCoeffs - 1);
  (*error).resize(rows, 1);
  
  (*lpcCoeffs).setZero();
  (*lpcCoeffs).col(0).setOnes();
  
  (*reflectionCoeffs).setZero();

  (*error).col(0) = _acorr.col(0);

  for ( int row = 0; row < rows; row++) {  
    Real gamma;
    
    for ( int i = 1; i < _numCoeffs; i++ ) {
      gamma = _acorr(row, i);
      
      // TODO: fix this when Eigen allows reverse()
      //if ( i >= 2) {
        //gamma += ((*lpcCoeffs).row(row).segment(1, i-1) * _acorr.row(row).segment(1, i-1).transpose().reverse())(0,0);
      //}
      
      for (int j = 1; j <= i-1; ++j) {
        gamma += (*lpcCoeffs)(row, j) * _acorr(row, i-j);  
      }
      
      (*reflectionCoeffs)(row, i-1) = - gamma / (*error)(row, 0);

      
      (*error)(row, 0) *= (1 - (*reflectionCoeffs)(row, i-1) * (*reflectionCoeffs).conjugate()(row, i-1));
      
      if(i >= 2){
        _temp = (*lpcCoeffs).block(row, 1, 1, i-1);
        reverseCols(&_temp);
        
        (*lpcCoeffs).block(row, 1, 1, i-1) += (*reflectionCoeffs)(row, i-1) * _temp.conjugate();
      }
      
      (*lpcCoeffs)(row, i) = (*reflectionCoeffs)(row, i-1);
    }
  }
  
  DEBUG("LPC: Finished Processing");
}

void LPC::reset(){
  // Initial values
}

int LPC::numCoeffs() const {
  return _numCoeffs;
}
