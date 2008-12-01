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
#include <vector>

#include "melbands.h"

#include "debug.h"
#include "typedefs.h"

using namespace std;

// define a custom template unary functor
template<typename Scalar>
struct CwiseCeilOp {
  CwiseCeilOp(){}
  const Scalar operator()(const Scalar& x) const { return ceil(x); }
};

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

MelBands::MelBands(Real lowFreq, Real highFreq, int numBands, Real samplerate, int spectrumLength) {
  DEBUG("MELBANDS: Constructor lowFreq: " << lowFreq << ", highFreq: " << highFreq << ", numBands: " << numBands << ", samplerate: " << samplerate << ", spectrumLength: " << spectrumLength);

  if ( lowFreq >= highFreq ) {
    // Throw an exception
  }

  if ( numBands <= 0 ) {
    // Throw an exception
  }
  
  _spectrumLength = spectrumLength;
  _samplerate = samplerate;
  _lowFreq = lowFreq;
  _highFreq = highFreq;
  _numBands = numBands;

  DEBUG("MELBANDS: Constructed");
}

void MelBands::setup(){
  DEBUG("MELBANDS: Setting up...");
  Real highMel = linearToMelRealStevens1937(_highFreq);
  Real lowMel = linearToMelRealStevens1937(_lowFreq);
  DEBUG("MELBANDS: lowMel: " << lowMel << ", highMel: " << highMel);

  Real stepMel = (highMel - lowMel) / (_numBands + 1.0);
  Real stepSpectrum = Real(_spectrumLength) / _samplerate;
  
  // start Mel frequencies of filters
  MatrixXR starts(_numBands, 1);
  for (int i=0; i<starts.rows(); i++) {
    starts(i, 0) = (Real(i) * stepMel + lowMel);
  }
  MatrixXR startsLinear = melToLinearStevens1937(starts) * stepSpectrum;

  // center Mel frequencies of filters
  MatrixXR centers(_numBands, 1);
  for (int i=0; i<centers.rows(); i++) {
    centers(i, 0) = (Real(i + 1) * stepMel + lowMel);
  }

  MatrixXR centersLinear = melToLinearStevens1937(centers) * stepSpectrum;

  // stop Mel frequencies of filters
  MatrixXR stops(_numBands, 1);
  for (int i=0; i<stops.rows(); i++) {
    stops(i, 0) = (Real(i + 2) * stepMel + lowMel);
  }

  MatrixXR stopsLinear = melToLinearStevens1937(stops) * stepSpectrum;
  
  // start bins of filters
  MatrixXi startBins = startsLinear.unaryExpr(CwiseCeilOp<Real>());

  // stop bins of filters
  MatrixXi stopBins = stopsLinear.unaryExpr(CwiseCeilOp<Real>());

  // set the start bins
  _starts.set(startBins);
  
  // fill in the weights
  for (int i=0; i < _starts.rows(); i++) {
    int startBin = startBins(i, 0);
    int stopBin = stopBins(i, 0);
    
    int filterLength = stopBin - startBin;
    
    DEBUG("MELBANDS: filterLength: " << filterLength);
    MatrixXR newFilter = MatrixXR::Constant(1, 1, 1.0);
    if (filterLength != 0){
      newFilter.resize(filterLength, 1);
      
      Real start = startsLinear(i, 0);
      Real center = centersLinear(i, 0);
      Real stop = stopsLinear(i, 0);
      
      triangleWindow(&newFilter, start - startBin, stop  - startBin, center  - startBin);
    }
    
    _weights.push_back(newFilter);
  }

  DEBUG("MELBANDS: Finished set up...");
}

void MelBands::triangleWindow(MatrixXR* window, Real start, Real stop, Real center, Real height) {
  int filterLength = (*window).rows();

  if (center == -1) {
    DEBUG("MELBANDS: Triangle window setting the center by default");
    center = start + (stop - start) / 2.0;
  }

  if ((center < 0) || (center > filterLength) || (center == start) || (center == stop)) {
    // Throw an exception invalid filter center
  }

  DEBUG("MELBANDS: Creating triangle window");
  DEBUG("MELBANDS: Triangle start: " << start << ", stop: " << stop << ", center: " << center << ", height: " << height);
  
  for (int i = 0; i < filterLength; i++) {
    if (i <= center) {
      (*window)(i,0) =  height * (Real(i) - start) / (center - start);
    } else {
      (*window)(i,0) =  height * (Real(1.0) - ((Real(i) - center) / (stop - center)));
    }
  }
  
  DEBUG("MELBANDS: Triangle window created: [" << (*window).transpose() << "]");
}

/**
 *
 * Mel scales computed using the original formula proposed by:
 *
 * Stevens, Stanley Smith; Volkman; John; & Newman, Edwin. (1937). 
 * A scale for the measurement of the psychological magnitude of pitch.
 * Journal of the Acoustical Society of America, 8 (3), 185-190.
 *
 */
Real MelBands::linearToMelRealStevens1937(Real linearFreq) {
  return log((linearFreq / 700.0) + 1.0) * 1127.01048;
}

Real MelBands::melToLinearRealStevens1937(Real melFreq) {
  return (exp(melFreq / 1127.01048) - 1.0) * 700.0;
}

MatrixXR MelBands::linearToMelStevens1937(MatrixXR linearFreq) {
  return ((linearFreq / 700.0).cwise() + 1.0).cwise().log() * 1127.01048;
}

MatrixXR MelBands::melToLinearStevens1937(MatrixXR melFreq) {
  return ((melFreq / 1127.01048).cwise().exp().cwise() - 1.0) * 700.0;
}


/**
 *
 * Mel scales computed using the formula proposed by:
 *  
 * Fant, Gunnar. (1968).
 * Analysis and synthesis of speech processes.
 * In B. Malmberg (Ed.), Manual of phonetics (pp. 173-177). Amsterdam: North-Holland.
 *
 */
Real MelBands::linearToMelRealFant1968(Real linearFreq) {
  return (1000.0 / log(2.0)) * log(1.0 + linearFreq/1000.0);
}

Real MelBands::melToLinearRealFant1968(Real melFreq) {
  return 1000.0 * (exp(melFreq * log(2.0) / 1000.0) - 1.0);
}

MatrixXR MelBands::linearToMelFant1968(MatrixXR linearFreq) {
  return (1000.0 / log(2.0)) * ((linearFreq / 1000.0).cwise() + 1.0).cwise().log();
}

MatrixXR MelBands::melToLinearFant1968(MatrixXR melFreq) {
  return 1000.0 * ((melFreq * log(2.0) / 1000.0).cwise().exp().cwise() - 1.0);
}

void MelBands::reset(){

}
