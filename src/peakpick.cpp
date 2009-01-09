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

#include <cmath>

#include "peakpick.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

PeakPick::PeakPick(int numPeaks) {
  DEBUG("PEAKPICK: Constructor numPeaks: " << numPeaks);
  
  _numPeaks = numPeaks;

  setup();
  DEBUG("PEAKPICK: Constructed");
}

PeakPick::~PeakPick() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void PeakPick::setup(){
  // Prepare the buffers
  DEBUG("PEAKPICK: Setting up...");

  reset();
  DEBUG("PEAKPICK: Finished set up...");
}


void PeakPick::process(MatrixXC fft, MatrixXR* peakPositions, MatrixXR* peakMagnitudes){
  DEBUG("PEAKPICK: Processing");
  int peakIndex;
  
  Real magnitude;
  Real pastMagnitude;
  Real postMagnitude;

  int numPeaks = _numPeaks;
  if(numPeaks == -1){
    numPeaks = fft.cols();
  }

  DEBUG("PEAKPICK: Processing, fft.shape: (" << fft.rows() << ", " << fft.cols() << ")");

  (*peakPositions).resize(fft.rows(), numPeaks);
  (*peakPositions).setConstant(-1);

  (*peakMagnitudes).resize(fft.rows(), numPeaks);
  (*peakMagnitudes).setConstant(-1);

  _magnitudes.set(fft.cwise().abs());

  DEBUG("PEAKPICK: Processing, _magnitudes.shape: (" << _magnitudes.rows() << ", " << _magnitudes.cols() << ")");


  for ( int j = 0 ; j < _magnitudes.rows(); j++){
    peakIndex = 0;
    
    magnitude = 0;
    pastMagnitude = 0;
    postMagnitude = 0;
  
    for ( int i = -1; i < _magnitudes.row(j).cols() - 1; i++) {
      if ( peakIndex >= _numPeaks ) {
        break;
      }

      pastMagnitude = magnitude;
      magnitude = postMagnitude;
      postMagnitude = _magnitudes( j, i + 1 );
      
      if ((magnitude > pastMagnitude) && (magnitude >= postMagnitude)) {
        (*peakMagnitudes)(j, peakIndex) = magnitude;
        (*peakPositions)(j, peakIndex) = i;
        peakIndex ++;
      }
    }
  }
  DEBUG("PEAKPICK: Finished Processing");
}

void PeakPick::reset(){
  // Initial values
}

int PeakPick::numPeaks() const {
  return _numPeaks;
}
