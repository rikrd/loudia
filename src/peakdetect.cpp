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

#include "peakdetect.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

PeakDetect::PeakDetect(int numPeaks, int minPeakWidth) {
  DEBUG("PEAKDETECT: Constructor numPeaks: " << numPeaks << ", minPeakWidth: " << minPeakWidth);
  
  _numPeaks = numPeaks;
  _minPeakWidth = minPeakWidth;

  setup();
  DEBUG("PEAKDETECT: Constructed");
}

PeakDetect::~PeakDetect() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void PeakDetect::setup(){
  // Prepare the buffers
  DEBUG("PEAKDETECT: Setting up...");

  reset();

  DEBUG("PEAKDETECT: Finished set up...");
}


void PeakDetect::process(const MatrixXC& fft, MatrixXR* peakPositions, MatrixXR* peakMagnitudes){
  DEBUG("PEAKDETECT: Processing");
  int peakIndex;
  
  int numPeaks = _numPeaks;
  if(numPeaks == -1){
    numPeaks = fft.cols();
  }

  DEBUG("PEAKDETECT: Processing, fft.shape: (" << fft.rows() << ", " << fft.cols() << ")");

  (*peakPositions).resize(fft.rows(), numPeaks);
  (*peakPositions).setConstant(-1);

  (*peakMagnitudes).resize(fft.rows(), numPeaks);
  (*peakMagnitudes).setConstant(-1);

  _magnitudes.set(fft.cwise().abs());

  DEBUG("PEAKDETECT: Processing, _magnitudes.shape: (" << _magnitudes.rows() << ", " << _magnitudes.cols() << ")");

  int maxRow;
  int maxCol;
  RowXR band(_minPeakWidth);
  
  for ( int i = 0 ; i < _magnitudes.rows(); i++){
    peakIndex = 0;
    
    band.setZero();

    for ( int j = (_minPeakWidth / 2); j < _magnitudes.row(i).cols() - (_minPeakWidth / 2); j++) {
      if ( peakIndex >= _numPeaks ) {
        break;
      }

      int inf = j - (_minPeakWidth / 2);

      band = _magnitudes.row(i).segment(inf, _minPeakWidth);

      band.maxCoeff( &maxRow, &maxCol );

      if ( maxCol == floor(_minPeakWidth / 2) ) {
        (*peakMagnitudes)(i, peakIndex) = _magnitudes(i, j);
        (*peakPositions)(i, peakIndex) = j;
        peakIndex ++;
      }
    }
  }
  DEBUG("PEAKDETECT: Finished Processing");
}

void PeakDetect::reset(){
  // Initial values
}

int PeakDetect::numPeaks() const {
  return _numPeaks;
}

int PeakDetect::minPeakWidth() const {
  return _minPeakWidth;
}
