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

  DEBUG("MELBANDS: Constructor {lowFreq:" << lowFreq << ", highFreq:" << highFreq << ", numBands:" << numBands << ", samplerate:" << samplerate << ", spectrumLength:" << spectrumLength << "}");

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
  
}

void MelBands::setup(){
  Real highMel = linearToMelReal(_highFreq);
  Real lowMel = linearToMelReal(_lowFreq);
  Real stepMel = (highMel - lowMel) / (_numBands + 1.0);
  Real stepSpectrum = Real(_spectrumLength) / _samplerate;
  
  // start Mel frequencies of filters
  MatrixXR starts(_numBands, 1);
  for (int i=0; i<starts.rows(); i++) {
    starts(i, 0) = (Real(i) * stepMel + lowMel);
  }
  MatrixXR startsLinear = melToLinear(starts) * stepSpectrum;

  // center Mel frequencies of filters
  MatrixXR centers(_numBands, 1);
  for (int i=0; i<centers.rows(); i++) {
    centers(i, 0) = (Real(i + 1) * stepMel + lowMel);
  }

  MatrixXR centersLinear = melToLinear(centers) * stepSpectrum;

  // stop Mel frequencies of filters
  MatrixXR stops(_numBands, 1);
  for (int i=0; i<stops.rows(); i++) {
    stops(i, 0) = (Real(i + 2) * stepMel + lowMel);
  }

  MatrixXR stopsLinear = melToLinear(stops) * stepSpectrum;
  
  // start bins of filters
  MatrixXi startBins = startsLinear.unaryExpr(CwiseCeilOp<Real>());

  // stop bins of filters
  MatrixXi stopBins = stopsLinear.unaryExpr(CwiseCeilOp<Real>());

  // set the start bins
  _starts = startBins;
  
  // fill in the weights
  for (int i=0; i < _starts.rows(); i++) {
    int startBin = startBins(i, 0);
    int stopBin = stopBins(i, 0);
    
    int filterLength = stopBin - startBin;

    MatrixXR newFilter(filterLength, 1);

    Real start = startsLinear(i, 0);
    Real center = centersLinear(i, 0);
    Real stop = stopsLinear(i, 0);
    
    triangleWindow(&newFilter, start - startBin, stop  - startBin, center  - startBin);
    _weights.push_back(newFilter);
  }
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

Real MelBands::linearToMelReal(Real linearFreq) {
  return log((linearFreq / 700) + 1.0) * 1127.0;
}

Real MelBands::melToLinearReal(Real melFreq) {
  return (exp(melFreq / 1127.0) - 1.0) * 700.0;
}

MatrixXR MelBands::linearToMel(MatrixXR linearFreq) {
  return ((linearFreq / 700).cwise() + 1.0).cwise().log() * 1127.0;
}

MatrixXR MelBands::melToLinear(MatrixXR melFreq) {
  return ((melFreq / 1127.0).cwise().exp().cwise() - 1.0) * 700.0;
}

void MelBands::reset(){

}
