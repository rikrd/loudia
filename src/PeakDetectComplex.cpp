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

struct byMagnitudeComp{
  bool operator() (peak i, peak j) { return ( i.mag > j.mag ); }
} byMagnitudeComplex;

struct byPositionComp{
  bool operator() (peak i, peak j) { return ( i.pos < j.pos ); }
} byPositionComplex;

PeakDetectComplex::PeakDetectComplex(int peakCount, SortMethod sortMethod, int minimumPeakWidth, int candidateCount, Real minimumPeakContrast) :
  _peakCount( peakCount ),
  _minimumPeakWidth( minimumPeakWidth ),
  _candidateCount( candidateCount ),
  _minimumPeakContrast( minimumPeakContrast ),
  _sortMethod( sortMethod )

{
  DEBUG("PEAKDETECTCOMPLEX: Constructor peakCount: " << _peakCount 
        << ", minimumPeakWidth: " << _minimumPeakWidth
        << ", candidateCount: " << _candidateCount);
  
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


void PeakDetectComplex::process(const MatrixXC& frames, 
                                MatrixXR* peakPositions, MatrixXR* peakMagnitudes, MatrixXR* peakPhases){
  DEBUG("PEAKDETECTCOMPLEX: Processing");
  
  const int rows = frames.rows();
  
  DEBUG("PEAKDETECTCOMPLEX: Processing, frames.shape: (" << rows << ", " << frames.cols() << ")");

  (*peakPositions).resize(rows, _peakCount);
  (*peakPositions).setConstant(-1);

  (*peakMagnitudes).resize(rows, _peakCount);
  (*peakMagnitudes).setConstant(-1);

  (*peakPhases).resize(rows, _peakCount);
  (*peakPhases).setConstant(-1);

  _magnitudes = frames.cwise().abs();
  _phases = frames.cwise().angle();

  DEBUG("PEAKDETECTCOMPLEX: Processing, _magnitudes.shape: (" << rows << ", " << _magnitudes.cols() << ")");
  
  int maxRow;
  int maxCol;
  
  Real maxVal;
  Real minVal;
  
  vector<peak> peaks;
  peaks.reserve(frames.cols());
  
  for ( int i = 0 ; i < rows; i++){

    peaks.clear();
    
    for ( int j = (_minimumPeakWidth / 2); j < _magnitudes.row(i).cols() - (_minimumPeakWidth / 2); j++) {
      // If we don't need sorting then only the first peakCount peaks are needed
      if( ( _sortMethod == NONE ) && ( (int)peaks.size() > _peakCount ) ) break;

      int inf = j - (_minimumPeakWidth / 2);
      
      // Get the maximum value and position of a region (corresponding to the min bandwidth of the peak)
      // of the spectrum
      maxVal = _magnitudes.row(i).segment(inf, _minimumPeakWidth).maxCoeff( &maxRow, &maxCol );
      
      // If the position of the maximum value is the center, then consider it as a peak candidate
      if ( maxCol == floor(_minimumPeakWidth / 2) ) {

        // Get the mininum value of the region
        minVal = _magnitudes.row(i).segment(inf, _minimumPeakWidth).minCoeff();

        // If the contrast is bigger than what minimumPeakContrast says, then select as peak
        if ( maxVal - minVal >= _minimumPeakContrast ) {

          peak p = {j, _magnitudes(i, j), _phases(i, j)};
          peaks.push_back(p);

        }
      }
    }
      
    // Get the largest candidates
    int candidateCount = (int)peaks.size();
    if(_candidateCount > 0) {
      candidateCount = min(candidateCount, _candidateCount);
      std::sort(peaks.begin(), peaks.end(), byMagnitudeComplex);
    }
    
    // Sort the candidates using position or magnitude
    switch ( _sortMethod ) {
    case BYPOSITION:      
      std::sort(peaks.begin(), peaks.begin() + candidateCount, byPositionComplex);
      break;
      
    case BYMAGNITUDE:
      // We have not done a candidate preselection, we must do the sorting
      if (_candidateCount <= 0)
        std::sort(peaks.begin(), peaks.begin() + candidateCount, byMagnitudeComplex);
      break;
      
    case NONE:
    default:
      break;
    }
    
    // Take the first peakCount
    int peakCount = min(_peakCount, candidateCount);      
    // Put the peaks in the matrices
    for( int j = 0; j < peakCount; j++ ){
      (*peakMagnitudes)(i, j) = peaks[j].mag;
      (*peakPositions)(i, j) = peaks[j].pos;
      (*peakPhases)(i, j) = peaks[j].phase;
    } 
  }

  DEBUG("PEAKDETECTCOMPLEX: Finished Processing");
}

void PeakDetectComplex::reset(){
  // Initial values
}

int PeakDetectComplex::peakCount() const {
  return _peakCount;
}

void PeakDetectComplex::setPeakCount( int count, bool callSetup ) {
  _peakCount = count;
  if ( callSetup ) setup();  
}

int PeakDetectComplex::candidateCount() const {
  return _candidateCount;
}

void PeakDetectComplex::setCandidateCount( int count, bool callSetup ) {
  _candidateCount = count;
  if ( callSetup ) setup();  
}

int PeakDetectComplex::minimumPeakWidth() const {
  return _minimumPeakWidth;
}

void PeakDetectComplex::setMinimumPeakWidth( int width, bool callSetup ) {
  _minimumPeakWidth = width;
  if ( callSetup ) setup();
}

int PeakDetectComplex::minimumPeakContrast() const {
  return _minimumPeakContrast;
}

void PeakDetectComplex::setMinimumPeakContrast( Real contrast, bool callSetup ) {
  _minimumPeakContrast = contrast;
  if ( callSetup ) setup();
}

PeakDetectComplex::SortMethod PeakDetectComplex::sortMethod() const {
  return _sortMethod;
}

void PeakDetectComplex::setSortMethod( SortMethod method, bool callSetup ) {
  _sortMethod = method;
  if ( callSetup ) setup();  
}
