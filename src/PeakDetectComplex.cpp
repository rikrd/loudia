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

#include "PeakDetectComplex.h"

using namespace std;
using namespace Eigen;

struct peak{
  Real pos;
  Real mag;
  Real phase;
  
  // A peak is smaller (first in the list)
  // if it's magnitude is larger
  bool operator <(peak const& other) const {
    return mag > other.mag;
  }
};

struct byMagnitudeComplexComp{
  bool operator() (peak i, peak j) { return ( i.mag > j.mag ); }
} byMagnitudeComplex;

struct byPositionComplexComp{
  bool operator() (peak i, peak j) { return ( i.pos < j.pos ); }
} byPositionComplex;

PeakDetectComplex::PeakDetectComplex(int numPeaks, SortType sort, int minPeakWidth, int numCandidates, Real minPeakContrast) :
  _numPeaks(numPeaks),
  _minPeakWidth(minPeakWidth),
  _numCandidates(numCandidates),
  _minPeakContrast(minPeakContrast),
  _sort(sort)

{
  DEBUG("PEAKDETECTCOMPLEX: Constructor numPeaks: " << _numPeaks 
        << ", minPeakWidth: " << _minPeakWidth
        << ", numCandidates: " << _numCandidates);
  
  setup();

  DEBUG("PEAKDETECTCOMPLEX: Constructed");
}

PeakDetectComplex::~PeakDetectComplex() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void PeakDetectComplex::setup(){
  // Prepare the buffers
  DEBUG("PEAKDETECTCOMPLEX: Setting up...");

  reset();

  DEBUG("PEAKDETECTCOMPLEX: Finished set up...");
}


void PeakDetectComplex::process(const MatrixXC& input, 
                                MatrixXR* peakPositions, MatrixXR* peakMagnitudes, MatrixXR* peakPhases){
  DEBUG("PEAKDETECTCOMPLEX: Processing");
  const int rows = input.rows();

  int numPeaks = _numPeaks;
  if( numPeaks == -1 ){
    numPeaks = input.cols();
  }
  
  DEBUG("PEAKDETECTCOMPLEX: Processing, input.shape: (" << input.rows() << ", " << input.cols() << ")");

  (*peakPositions).resize(rows, numPeaks);
  (*peakPositions).setConstant(-1);

  (*peakMagnitudes).resize(rows, numPeaks);
  (*peakMagnitudes).setConstant(-1);

  (*peakPhases).resize(rows, numPeaks);
  (*peakPhases).setConstant(-1);

  _magnitudes = input.cwise().abs();
  _phases = input.cwise().angle();

  DEBUG("PEAKDETECTCOMPLEX: Processing, _magnitudes.shape: (" << _magnitudes.rows() << ", " << _magnitudes.cols() << ")");
  
  int maxRow;
  int maxCol;
  
  Real maxVal;
  Real minVal;

  vector<peak> peaks;
  peaks.reserve(input.cols());

  for ( int i = 0 ; i < rows; i++){
    DEBUG("PEAKDETECTCOMPLEX: Processing, new row");
    peaks.clear();
    
    for ( int j = (_minPeakWidth / 2); j < _magnitudes.row(i).cols() - (_minPeakWidth / 2); j++) {
      // If we don't need sorting then only the first numPeaks peaks are needed
      if( ( _sort == NOSORT ) && ( (int)peaks.size() > numPeaks ) ) break;

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
        }
      }
    }

    DEBUG("PEAKDETECTCOMPLEX: Processing, get largest candidates");
    // Get the largest candidates
    int candidateCount = peaks.size();
    if(_numCandidates > 0) {
      candidateCount = min(candidateCount, _numCandidates);
      std::sort(peaks.begin(), peaks.end(), byMagnitudeComplex);
    }
    
    // Sort the candidates using position or magnitude
    switch ( _sort ) {
    case BYPOSITION:      
      std::sort(peaks.begin(), peaks.begin() + candidateCount, byPositionComplex);
      break;
      
    case BYMAGNITUDE:
      // We have not done a candidate preselection, we must do the sorting
      if (_numCandidates <= 0)
        std::sort(peaks.begin(), peaks.begin() + candidateCount, byMagnitudeComplex);
      break;
      
    case NOSORT:
    default:
      break;
    }
    
    DEBUG("PEAKDETECTCOMPLEX: Processing, take first numPeaks");      
    // Take the first numPeaks
    int peakCount = min(numPeaks, candidateCount);
    // Put the peaks in the matrices
    for( int j = 0; j < peakCount; j++ ){
      (*peakMagnitudes)(i, j) = peaks[j].mag;
      (*peakPhases)(i, j) = peaks[j].phase;
      (*peakPositions)(i, j) = peaks[j].pos;
    }

    DEBUG("PEAKDETECTCOMPLEX: Processing, finished row");
  }

  DEBUG("PEAKDETECTCOMPLEX: Finished Processing");
}

void PeakDetectComplex::reset(){
  // Initial values
}

int PeakDetectComplex::numPeaks() const {
  return _numPeaks;
}

int PeakDetectComplex::minPeakWidth() const {
  return _minPeakWidth;
}
