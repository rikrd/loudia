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

#include "PeakInterpolation.h"

using namespace std;
using namespace Eigen;

PeakInterpolation::PeakInterpolation() {
  LOUDIA_DEBUG("PEAKINTERPOLATION: Constructor");
  
  setup();

  LOUDIA_DEBUG("PEAKINTERPOLATION: Constructed");
}

PeakInterpolation::~PeakInterpolation() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void PeakInterpolation::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("PEAKINTERPOLATION: Setting up...");

  reset();

  LOUDIA_DEBUG("PEAKINTERPOLATION: Finished set up...");
}


void PeakInterpolation::process(const MatrixXR& input,
                              const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes,
                              MatrixXR* peakPositionsInterp, MatrixXR* peakMagnitudesInterp) {
  
  LOUDIA_DEBUG("PEAKINTERPOLATION: Processing");  
  Real leftMag;
  Real rightMag;
  Real mag, interpFactor;
  
  (*peakPositionsInterp).resize(input.rows(), peakPositions.cols());
  (*peakMagnitudesInterp).resize(input.rows(), peakPositions.cols());
  
  _magnitudes = input.cwise().abs();
  
  for ( int row = 0 ; row < _magnitudes.rows(); row++ ) {
  
    for ( int i = 0; i < peakPositions.cols(); i++ ) {
      
      // If the position is -1 do nothing since it means it is nothing
      if( peakPositions(row, i) == -1 ){

        (*peakMagnitudesInterp)(row, i) = peakMagnitudes(row, i);
        (*peakPositionsInterp)(row, i) = peakPositions(row, i);
        
      } else {
        
        // Take the center magnitude in dB
        mag = 20.0 * log10( peakMagnitudes(row, i) );

        // Take the left magnitude in dB
        if( peakPositions(row, i) <= 0 ){
          
          leftMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) + 1) );
          
        } else {
          
          leftMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) - 1) );
          
        }
        
        // Take the right magnitude in dB
        if( peakPositions(row, i) >= _magnitudes.row(row).cols() - 1 ){
          
          rightMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) - 1) );
          
        } else {
          
          rightMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) + 1) );
          
        }
                
        // Calculate the interpolated position
        (*peakPositionsInterp)(row, i) = peakPositions(row, i) + 0.5 * (leftMag - rightMag) / (leftMag - 2.0 * mag + rightMag);

        interpFactor = ((*peakPositionsInterp)(row, i) - peakPositions(row, i));

        // Calculate the interpolated magnitude in dB
        (*peakMagnitudesInterp)(row, i) = mag - 0.25 * (leftMag - rightMag) * interpFactor;
      }
    }
  }
  
  LOUDIA_DEBUG("PEAKINTERPOLATION: Finished Processing");
}

void PeakInterpolation::reset(){
  // Initial values
}

