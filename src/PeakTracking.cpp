/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
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

#include "Typedefs.h"
#include "Debug.h"

#include "PeakTracking.h"
#include <limits>

using namespace std;
using namespace Eigen;

PeakTracking::PeakTracking(int trajectoryCount, Real maximumFrequencyChange, int silentFrameCount)
{
  DEBUG("PEAKTRACKING: Constructor");
  
  setTrajectoryCount( trajectoryCount, false );
  setMaximumFrequencyChange( maximumFrequencyChange, false );
  setSilentFrameCount( silentFrameCount, false );

  setup();
  
  DEBUG("PEAKTRACKING: Constructed");
}

PeakTracking::~PeakTracking() {
}


void PeakTracking::setup(){
  // Prepare the buffers
  DEBUG("PEAKTRACKING: Setting up...");
  
  reset();
  
  DEBUG("PEAKTRACKING: Finished set up...");
}


void PeakTracking::process(const MatrixXC& fft,
                           const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes,
                           MatrixXR* trajPositions, MatrixXR* trajMagnitudes){
  
  DEBUG("PEAKTRACKING: Processing");  
  
  (*trajPositions).resize(fft.rows(), _trajectoryCount);
  (*trajMagnitudes).resize(fft.rows(), _trajectoryCount);
  
  (*trajPositions) = MatrixXR::Constant(fft.rows(), _trajectoryCount, -1.0);
  (*trajMagnitudes) = MatrixXR::Constant(fft.rows(), _trajectoryCount, -120.0);

  MatrixXR currPeakPositions = peakPositions;
  MatrixXR currPeakMagnitudes = peakMagnitudes;
  
  for ( int row = 0 ; row < fft.rows(); row++ ) {

    // Find the closest peak to each of the trajectories
    for ( int i = 0 ; i < _pastTrajPositions.cols(); i++  ) {
      
      if( ! isinf( _pastTrajPositions(row, i) ) ) {
        
        int posRow, posCol;
        Real minFreqBinChange = (currPeakPositions.row(row).cwise() - _pastTrajPositions(row, i)).cwise().abs().minCoeff(&posRow, &posCol);
        
        if ( minFreqBinChange <= _maximumFrequencyChange ) {
          // A matching peak has been found
          DEBUG("PEAKTRACKING: Processing 'Matching peak: " << posCol << "' minFreqBinChange: " << minFreqBinChange);
          
          (*trajPositions)(row, i) = currPeakPositions(row, posCol);
          (*trajMagnitudes)(row, i) = currPeakMagnitudes(row, posCol);
          
          _pastTrajPositions(0, i) = (*trajPositions)(row, i);
          _pastTrajMagnitudes(0, i) = (*trajMagnitudes)(row, i);

          currPeakPositions(row, posCol) = numeric_limits<Real>::infinity();
          currPeakMagnitudes(row, posCol) = numeric_limits<Real>::infinity();
          
        } else {
          // No matching peak has been found
          DEBUG("PEAKTRACKING: Processing 'No matching peaks' minFreqBinChange: " << minFreqBinChange);
          
          if ( _pastTrajMagnitudes(0, i) <= (-120.0 - _silentFrameCount) ) {

            // The trajectory has been silent too long (resetting it)

            _pastTrajMagnitudes(0, i) = numeric_limits<Real>::infinity();
            _pastTrajPositions(0, i) = numeric_limits<Real>::infinity();

          } else if ( _pastTrajMagnitudes(0, i) <= -120.0 ) {

            // The trajectory has been silent for one more frame

            _pastTrajMagnitudes(0, i) -= 1;

          } else {

            // The first frame the trajectory is silent

            _pastTrajMagnitudes(0, i) = -120.0;
            
          }
          
          (*trajPositions)(row, i) = isinf(_pastTrajPositions(0, i)) ? -1 : _pastTrajPositions(0, i);
          (*trajMagnitudes)(row, i) = -120.0;
          
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
          DEBUG("PEAKTRACKING: Processing the trajectory could not be created");
        }
      }  
    }
  }
  
  DEBUG("PEAKTRACKING: Finished Processing");
}

bool PeakTracking::createTrajectory(Real peakPos, Real peakMag,
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

void PeakTracking::reset(){
  // Initial values
  if ( !numeric_limits<Real>::has_infinity ) {
    // Throw PlatformError infinity not supported
  }

  Real inf = numeric_limits<Real>::infinity();
  _pastTrajPositions = MatrixXR::Constant(1, _trajectoryCount, inf);
  _pastTrajMagnitudes = MatrixXR::Constant(1, _trajectoryCount, inf);
}

int PeakTracking::trajectoryCount() const {
  return _trajectoryCount;
}

void PeakTracking::setTrajectoryCount( int count, bool callSetup ) {
  _trajectoryCount = count;
  if ( callSetup ) setup();  
}

Real PeakTracking::maximumFrequencyChange() const {
  return _maximumFrequencyChange;
}

void PeakTracking::setMaximumFrequencyChange( Real change, bool callSetup ) {
  _maximumFrequencyChange = change;
  if ( callSetup ) setup();
}

int PeakTracking::silentFrameCount() const {
  return _silentFrameCount;
}

void PeakTracking::setSilentFrameCount( int count, bool callSetup ) {
  _silentFrameCount = count;
  if ( callSetup ) setup();  
}
