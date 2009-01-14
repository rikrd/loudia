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

#include "peakcontinue.h"
#include <cmath>
#include <limits>

using namespace std;

// import most common Eigen types 
using namespace Eigen;

PeakContinue::PeakContinue(int numTrajectories, Real maxFreqBinChange) {
  DEBUG("PEAKCONTINUE: Constructor");
  
  _numTrajectories = numTrajectories;
  _maxFreqBinChange = maxFreqBinChange;
  
  setup();
  
  DEBUG("PEAKCONTINUE: Constructed");
}

PeakContinue::~PeakContinue() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void PeakContinue::setup(){
  // Prepare the buffers
  DEBUG("PEAKCONTINUE: Setting up...");
  
  reset();
  
  DEBUG("PEAKCONTINUE: Finished set up...");
}


void PeakContinue::process(const MatrixXC& fft,
                           const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes,
                           MatrixXR* trajPositions, MatrixXR* trajMagnitudes){
  
  DEBUG("PEAKCONTINUE: Processing");  
  
  (*trajPositions).resize(fft.rows(), _numTrajectories);
  (*trajMagnitudes).resize(fft.rows(), _numTrajectories);
  
  (*trajPositions).set(MatrixXR::Constant(fft.rows(), _numTrajectories, -1.0));
  (*trajMagnitudes).set(MatrixXR::Constant(fft.rows(), _numTrajectories, -120.0));

  MatrixXR currPeakPositions = peakPositions;
  MatrixXR currPeakMagnitudes = peakMagnitudes;

  for ( int row = 0 ; row < fft.rows(); row++ ) {

    // Find the closest peak to each of the trajectories
    for ( int i = 0 ; i < _pastTrajPositions.cols(); i++  ) {
      
      if( ! isinf( _pastTrajPositions(row, i) ) ) {
        
        int posRow, posCol;
        Real minFreqBinChange = (currPeakPositions.row(row).cwise() - _pastTrajPositions(row, i)).cwise().abs().minCoeff(&posRow, &posCol);
        
        if (minFreqBinChange <= _maxFreqBinChange) {
          // A matching peak has been found
          DEBUG("PEAKCONTINUE: Processing 'Matching peak: " << posCol << "' minFreqBinChange: " << minFreqBinChange);
          
          if ( isinf( currPeakMagnitudes(row, posCol) ) ) {
            cout << "ERROR !!!!!!!!!!!!" << endl;
          }

          (*trajPositions)(row, i) = currPeakPositions(row, posCol);
          (*trajMagnitudes)(row, i) = currPeakMagnitudes(row, posCol);
          
          _pastTrajPositions(0, i) = (*trajPositions)(row, i);
          _pastTrajMagnitudes(0, i) = (*trajMagnitudes)(row, i);

          currPeakPositions(row, posCol) = numeric_limits<Real>::infinity();
          currPeakMagnitudes(row, posCol) = numeric_limits<Real>::infinity();
          
        } else {
          // No matching peak has been found
          DEBUG("PEAKCONTINUE: Processing 'No matching peaks' minFreqBinChange: " << minFreqBinChange);
          
          //_pastTrajPositions(0, i) = numeric_limits<Real>::infinity();
          //_pastTrajMagnitudes(0, i) = numeric_limits<Real>::infinity();
          (*trajPositions)(row, i) = _pastTrajPositions(0, i);
          (*trajMagnitudes)(row, i) = -120.0;
          
          _pastTrajMagnitudes(0, i) = -120.0;
          
        }
      }
    }
      
    // Find those peaks that haven't been assigned and create new trajectories
    for ( int i = 0; i < currPeakPositions.cols(); i++ ) {    
      Real pos = currPeakPositions(row, i);
      Real mag = currPeakMagnitudes(row, i);
        
      if( ! isinf( pos ) ){
        bool created = createTrajectory(pos, mag, 
                                        &_pastTrajPositions, &_pastTrajMagnitudes,
                                        trajPositions, trajMagnitudes,
                                        row);
        
        if (! created ){
          DEBUG("PEAKCONTINUE: Processing the trajectory could not be created");
        }
      }  
    }
  }
  
  DEBUG("PEAKCONTINUE: Finished Processing");
}

bool PeakContinue::createTrajectory(Real peakPos, Real peakMag,
                                    MatrixXR* pastTrajPositions, MatrixXR* pastTrajMagnitudes,
                                    MatrixXR* trajPositions, MatrixXR* trajMagnitudes,
                                    int row) {

  int maxRow, maxCol;
  Real maxPos = (*pastTrajPositions).row(row).maxCoeff(&maxRow, &maxCol);

  if ( isinf( maxPos ) ) {
    //DEBUG("MAXCOL: " << maxCol);
    //DEBUG("ROW: " << row);
    
    (*pastTrajPositions)(0, maxCol) = peakPos;
    (*pastTrajMagnitudes)(0, maxCol) = peakMag;

    //DEBUG("Past: ");

    (*trajPositions)(row, maxCol) = peakPos;
    (*trajMagnitudes)(row, maxCol) = peakMag;

    return true;
  }

  return false;
}

void PeakContinue::reset(){
  // Initial values
  if ( !numeric_limits<Real>::has_infinity ) {
    // Throw PlatformError infinity not supported
  }

  Real inf = numeric_limits<Real>::infinity();
  _pastTrajPositions.set(MatrixXR::Constant(1, _numTrajectories, inf));
  _pastTrajMagnitudes.set(MatrixXR::Constant(1, _numTrajectories, inf));
}

