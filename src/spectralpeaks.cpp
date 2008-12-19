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

#include "spectralpeaks.h"

#include "typedefs.h"
#include "debug.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

SpectralPeaks::SpectralPeaks(int numPeaks) {
  DEBUG("SpectralPeaks: Constructor lowFreq: " << lowFreq << ", highFreq: " << highFreq << ", numBands: " << numBands << ", samplerate: "<< samplerate << ", spectrumLength: " << spectrumLength << ", numCoeffs: " << numCoeffs);
  
  _numPeaks = numPeaks;

}

SpectralPeaks::~SpectralPeaks() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void SpectralPeaks::setup(){
  // Prepare the buffers
  DEBUG("SpectralPeaks: Setting up...");

  reset();
  DEBUG("SpectralPeaks: Finished set up...");
}


void SpectralPeaks::process(MatrixXR spectrum, MatrixXR* peakMagnitudes, MatrixXR* peakPositions){
  DEBUG("SpectralPeaks: Processing");
  int peakIndex;
  
  Real magnitude;
  Real pastMagnitude;
  Real postMagnitude;

  for ( int j = 0 ; j < spectrum.rows(); j++){
    peakIndex = 0;
    
    magnitude = 0;
    pastMagnitude = 0;
    postMagnitude = 0;
  
    for ( int i = 0; i < spectrum.row(j).cols() - 1; i++) {
      pastMagnitude = magnitude;
      magnitude = posMagnitude;
      postMagnitude = spectrum( j, i + 1 );
      
      if ((magnitude >= pastMagnitude) && (magnitude > postMagnitude)) {
        (*peakMagnitudes)(j, peakIndex) = magnitude;
        (*peakPositions)(j, peakIndex) = i;
      }
    }
  }
  DEBUG("SpectralPeaks: Finished Processing");
}

void SpectralPeaks::reset(){
  // Initial values
}

int SpectralPeaks::numPeaks() const {
  return _numPeaks;
}
