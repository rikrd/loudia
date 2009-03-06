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

#include "Typedefs.h"
#include "Debug.h"

#include <algorithm>
#include <vector>

#include "PeakDetect.h"

using namespace std;
using namespace Eigen;

struct peak{
  Real pos;
  Real mag;
  
  // A peak is smaller (first in the list)
  // if it's magnitude is larger
  bool operator <(peak const& other) const {
    return mag > other.mag;
  }
};

struct byMagnitudeComp{
  bool operator() (peak i, peak j) { return ( i.mag > j.mag ); }
} byMagnitude;

struct byPositionComp{
  bool operator() (peak i, peak j) { return ( i.pos < j.pos ); }
} byPosition;

PeakDetect::PeakDetect(int numPeaks, SortType sort, int minPeakWidth, int numCandidates, Real minPeakContrast) :
  _numPeaks(numPeaks),
  _minPeakWidth(minPeakWidth),
  _numCandidates(numCandidates),
  _minPeakContrast(minPeakContrast),
  _sort(sort)

{
  DEBUG("PEAKDETECT: Constructor numPeaks: " << numPeaks << ", minPeakWidth: " << minPeakWidth);
  
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
                         MatrixXR* peakPositions, MatrixXR* peakMagnitudes){
  DEBUG("PEAKDETECT: Processing");
  
  int numPeaks = _numPeaks;
  if( numPeaks == -1 ){
    numPeaks = fft.cols();
  }
  
  DEBUG("PEAKDETECT: Processing, fft.shape: (" << fft.rows() << ", " << fft.cols() << ")");

  (*peakPositions).resize(fft.rows(), numPeaks);
  (*peakPositions).setConstant(-1);

  (*peakMagnitudes).resize(fft.rows(), numPeaks);
  (*peakMagnitudes).setConstant(-1);

  _magnitudes = fft.cwise().abs();

  DEBUG("PEAKDETECT: Processing, _magnitudes.shape: (" << _magnitudes.rows() << ", " << _magnitudes.cols() << ")");
  
  int maxRow;
  int maxCol;
  
  Real maxVal;
  Real minVal;
  
  int peakIndex;
  
  for ( int i = 0 ; i < _magnitudes.rows(); i++){
    vector<peak> peaks;

    peakIndex = 0;

    // If we don't need sorting then only the first numPeaks peaks are needed
    if( ( _sort == NOSORT ) && ( peakIndex > numPeaks ) ) break;
    
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

          peak p = {j, _magnitudes(i, j)};
          peaks.push_back(p);

          peakIndex ++;
        }
      }      
      
      // Get the largest candidates
      int candidateCount = peakIndex;
      if(_numCandidates > 0) {
        candidateCount = min(peakIndex, _numCandidates);
        std::sort(peaks.begin(), peaks.end(), byMagnitude);
      }

      // Sort the candidates using position or magnitude
      switch ( _sort ) {
      case BYPOSITION:      
        std::sort(peaks.begin(), peaks.begin() + candidateCount, byPosition);
        break;

      case BYMAGNITUDE:
        if (_numCandidates <= 0)
          std::sort(peaks.begin(), peaks.begin() + candidateCount, byMagnitude);
        break;
        
      case NOSORT:
      default:
        break;
      }
      
      // Take the first numPeaks
      int peakCount = min(_numPeaks, candidateCount);      
      // Put the peaks in the matrices
      for( int j = 0; j < peakCount; j++ ){
        (*peakMagnitudes)(i, j) = peaks[j].mag;
        (*peakPositions)(i, j) = peaks[j].pos;
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
