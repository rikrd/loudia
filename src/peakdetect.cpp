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

#include <algorithm>
#include <vector>

#include "peakdetect.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

struct peak{
  Real pos;
  Real mag;
  Real phase;
  
  bool operator <(peak const& other) const {
    return mag > other.mag;
  }
};


PeakDetect::PeakDetect(int numPeaks, int minPeakWidth, Real minPeakContrast) {
  DEBUG("PEAKDETECT: Constructor numPeaks: " << numPeaks << ", minPeakWidth: " << minPeakWidth);
  
  _numPeaks = numPeaks;
  _minPeakWidth = minPeakWidth;
  _minPeakContrast = minPeakContrast;

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


void PeakDetect::process(const MatrixXC& fft, 
                         MatrixXR* peakPositions, MatrixXR* peakMagnitudes, MatrixXR* peakPhases){
  DEBUG("PEAKDETECT: Processing");
  
  int numPeaks = _numPeaks;
  if(numPeaks == -1){
    numPeaks = fft.cols();
  }
  
  DEBUG("PEAKDETECT: Processing, fft.shape: (" << fft.rows() << ", " << fft.cols() << ")");

  (*peakPositions).resize(fft.rows(), numPeaks);
  (*peakPositions).setConstant(-1);

  (*peakMagnitudes).resize(fft.rows(), numPeaks);
  (*peakMagnitudes).setConstant(-1);

  (*peakPhases).resize(fft.rows(), numPeaks);
  (*peakPhases).setConstant(-1);

  _magnitudes = fft.cwise().abs();
  _phases = fft.cwise().angle();

  DEBUG("PEAKDETECT: Processing, _magnitudes.shape: (" << _magnitudes.rows() << ", " << _magnitudes.cols() << ")");
  
  int maxRow;
  int maxCol;
  
  Real maxVal;
  Real minVal;
  
  int peakIndex;
  
  vector<peak> peaks;
  for ( int i = 0 ; i < _magnitudes.rows(); i++){
    peakIndex = 0;
    
    for ( int j = (_minPeakWidth / 2); j < _magnitudes.row(i).cols() - (_minPeakWidth / 2); j++) {     
      int inf = j - (_minPeakWidth / 2);
      
      // Get the maximum value and position of a region (corresponding to the min bandwidth of the peak)
      // of the spectrum
      maxVal = _magnitudes.row(i).segment(inf, _minPeakWidth).maxCoeff( &maxRow, &maxCol );
      
      // If the position of the maximum value is the center, then consider it as a peak candidate
      if ( maxCol == floor(_minPeakWidth / 2) ) {

        // Get the mininum value of the region
        minVal = _magnitudes.row(i).segment(inf, _minPeakWidth).minCoeff();

        // If the contrast is bigger than what minPeakContrast says, then select as peak
        if ( maxVal - minVal >= _minPeakContrast ) {

          peak p = {j, _magnitudes(i, j), _phases(i, j)};
          peaks.push_back(p);

          peakIndex ++;
        }
      }
    }
    
    // Order the peaks by magnitude
    std::sort(peaks.begin(), peaks.end());  

    int peakCount = min(_numPeaks, peakIndex);

    // Put the peaks in the matrices
    for( int j = 0; j < peakCount; j++ ){
      (*peakMagnitudes)(i, j) = peaks[j].mag;
      (*peakPhases)(i, j) = peaks[j].phase;
      (*peakPositions)(i, j) = peaks[j].pos;
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
