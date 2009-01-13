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
#include <vector>

#include "melbands.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

MelBands::MelBands(Real lowFreq, Real highFreq, int numBands, Real samplerate, int fftLength) 
{
  DEBUG("MELBANDS: Constructor lowFreq: " << lowFreq << ", highFreq: " << highFreq << ", numBands: " << numBands << ", samplerate: " << samplerate << ", fftLength: " << fftLength);

  if ( lowFreq >= highFreq ) {
    // Throw an exception
  }

  if ( numBands <= 0 ) {
    // Throw an exception
  }
  
  _fftLength = fftLength;
  _samplerate = samplerate;
  _lowFreq = lowFreq;
  _highFreq = highFreq;
  _numBands = numBands;

  setup();
  DEBUG("MELBANDS: Constructed");
}

void MelBands::setup(){
  DEBUG("MELBANDS: Setting up...");
    
  Real highMel = linearToMelGreenwood1990(_highFreq);
  Real lowMel = linearToMelGreenwood1990(_lowFreq);
  
  DEBUG("MELBANDS: lowMel: " << lowMel << ", highMel: " << highMel);

  Real stepMel = (highMel - lowMel) / (_numBands + 1.0);
  Real stepSpectrum = Real(_fftLength) / _samplerate;
  
  // start Mel frequencies of filters
  MatrixXR starts(_numBands, 1);
  for (int i=0; i<starts.rows(); i++) {
    starts(i, 0) = (Real(i) * stepMel + lowMel);
  }

  MatrixXR startsLinear;
  melToLinearMatrixGreenwood1990(starts, &startsLinear);
  startsLinear *= stepSpectrum;

  // stop Mel frequencies of filters
  MatrixXR stops(_numBands, 1);
  for (int i=0; i<stops.rows(); i++) {
    stops(i, 0) = (Real(i + 2) * stepMel + lowMel);
  }

  MatrixXR stopsLinear;
  melToLinearMatrixGreenwood1990(stops, &stopsLinear);
  stopsLinear *= stepSpectrum;


  // center Mel frequencies of filters
  MatrixXR centers(_numBands, 1);
  for (int i=0; i<centers.rows(); i++) {
    centers(i, 0) = (Real(i + 1) * stepMel + lowMel);
  }

  MatrixXR centersLinear = startsLinear + (stopsLinear - startsLinear) / 2.0;
  //melToLinearMatrixGreenwood1990(centers, &centersLinear);
  //centersLinear *= stepSpectrum;
  
  // start bins of filters
  MatrixXi startBins = startsLinear.cwise().ceil();

  // stop bins of filters
  MatrixXi stopBins = stopsLinear.cwise().ceil();

  std::vector<MatrixXR> weights;

  // fill in the weights
  for (int i=0; i < startBins.rows(); i++) {
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
    
    weights.push_back(newFilter);
  }

  _bands.setStartsWeights(startBins, weights);

  DEBUG("MELBANDS: Finished set up...");
}

void MelBands::process(MatrixXR spectrum, MatrixXR* bands) {  
  _bands.process(spectrum, bands);
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
 * Mel scales computed using the Greenwood function:
 *
 * Greenwood, DD. (1990)
 * A cochlear frequency-position function for several species - 29 years later,
 * Journal of the Acoustical Society of America, vol. 87, pp. 2592-2605.
 *
 */
Real MelBands::linearToMelGreenwood1990(Real linearFreq) {
  return log10((linearFreq / 165.4) + 1.0) / 2.1;
}

Real MelBands::melToLinearGreenwood1990(Real melFreq) {
  return 165.4 * (pow(10.0, 2.1 * melFreq) - 1.0);
}

void MelBands::linearToMelMatrixGreenwood1990(MatrixXR linearFreq, MatrixXR* melFreq) {
  DEBUG("MELBANDS: Scaling (Greenwood 1990) linearFreq: " << linearFreq);

  (*melFreq).set(((linearFreq / 165.4).cwise() + 1.0).cwise().logN(10) / 2.1);
}

void MelBands::melToLinearMatrixGreenwood1990(MatrixXR melFreq, MatrixXR* linearFreq) {
  DEBUG("MELBANDS: Scaling (Greenwood 1990) melFreq: " << melFreq);

  (*linearFreq).set(165.4 * ((melFreq * 2.1).cwise().expN(10.0).cwise() - 1.0));
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
Real MelBands::linearToMelStevens1937(Real linearFreq) {
  return log((linearFreq / 700.0) + 1.0) * 1127.01048;
}

Real MelBands::melToLinearStevens1937(Real melFreq) {
  return (exp(melFreq / 1127.01048) - 1.0) * 700.0;
}

void MelBands::linearToMelMatrixStevens1937(MatrixXR linearFreq, MatrixXR* melFreq) {
  (*melFreq).set(((linearFreq / 700.0).cwise() + 1.0).cwise().log() * 1127.01048);
}

void MelBands::melToLinearMatrixStevens1937(MatrixXR melFreq, MatrixXR* linearFreq) {
  (*linearFreq).set(((melFreq / 1127.01048).cwise().exp().cwise() - 1.0) * 700.0);
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
Real MelBands::linearToMelFant1968(Real linearFreq) {
  return (1000.0 / log(2.0)) * log(1.0 + linearFreq / 1000.0);
}

Real MelBands::melToLinearFant1968(Real melFreq) {
  return 1000.0 * (exp(melFreq * log(2.0) / 1000.0) - 1.0);
}

void MelBands::linearToMelMatrixFant1968(MatrixXR linearFreq, MatrixXR* melFreq) {
  (*melFreq).set((1000.0 / log(2.0)) * ((linearFreq / 1000.0).cwise() + 1.0).cwise().log());
}

void MelBands::melToLinearMatrixFant1968(MatrixXR melFreq, MatrixXR* linearFreq) {
  (*linearFreq).set(1000.0 * ((melFreq * log(2.0) / 1000.0).cwise().exp().cwise() - 1.0));
}

void MelBands::reset(){
  // Initial values
  _bands.reset();
}

void MelBands::starts(MatrixXI* result) const {
  return _bands.starts( result );
}

std::vector<MatrixXR> MelBands::weights() const {
  return _bands.weights();
}

void MelBands::bandWeights(int band, MatrixXR* bandWeights) const {
  return _bands.bandWeights( band, bandWeights );
}

int MelBands::bands() const {
  return _bands.bands();
}
