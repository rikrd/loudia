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

#include "MelBands.h"
#include "MelScales.h"

#include <vector>

using namespace std;
using namespace Eigen;

MelBands::MelBands(Real lowFrequency, Real highFrequency, int bandCount, Real samplerate, int fftSize, ScaleType scaleType) 
{
  
  DEBUG("MELBANDS: Constructor lowFrequency: " << _lowFrequency << 
        ", highFrequency: " << _highFrequency << 
        ", bandCount: " << _bandCount << 
        ", samplerate: " << _samplerate << 
        ", fftSize: " << _fftSize << 
        ", scaleType:" << _scaleType);

  if ( lowFrequency >= highFrequency ) {
    // Throw an exception, highFrequency must be higher than lowFrequency
  }

  if ( bandCount <= 0 ) {
    // Throw an exception, bandCount must be higher than 0
  }
  
  setLowFrequency( lowFrequency, false );
  setHighFrequency( highFrequency, false );
  setBandCount( bandCount, false );
  setSamplerate( samplerate, false );
  setFftSize( fftSize, false );
  setScaleType( scaleType, false );
  
  setup();
  
  DEBUG("MELBANDS: Constructed");
  
}

void MelBands::setup(){
  DEBUG("MELBANDS: Setting up...");
  
  // Set the linearToMel and melToLinear functions
  switch(_scaleType) {
  case STEVENS:
    _linearToMel = &linearToMelStevens1937;
    _melToLinear = &melToLinearStevens1937;
  
    _linearToMelMatrix = &linearToMelMatrixStevens1937;
    _melToLinearMatrix = &melToLinearMatrixStevens1937;

    break;

  case FANT:
    _linearToMel = &linearToMelFant1968;
    _melToLinear = &melToLinearFant1968;
  
    _linearToMelMatrix = &linearToMelMatrixFant1968;
    _melToLinearMatrix = &melToLinearMatrixFant1968;

    break;

  case GREENWOOD:
    _linearToMel = &linearToMelGreenwood1990;
    _melToLinear = &melToLinearGreenwood1990;
  
    _linearToMelMatrix = &linearToMelMatrixGreenwood1990;
    _melToLinearMatrix = &melToLinearMatrixGreenwood1990;

    break;

  }
  
  Real highMel = _linearToMel( _highFrequency );
  Real lowMel = _linearToMel( _lowFrequency );
  
  DEBUG("MELBANDS: lowMel: " << lowMel << ", highMel: " << highMel);

  Real stepMel = (highMel - lowMel) / (_bandCount + 1.0);
  Real stepSpectrum = Real(_fftSize) / _samplerate;
  
  // start Mel frequencies of filters
  MatrixXR starts(_bandCount, 1);
  for (int i=0; i<starts.rows(); i++) {
    starts(i, 0) = (Real(i) * stepMel + lowMel);
  }

  MatrixXR startsLinear;
  _melToLinearMatrix(starts, &startsLinear);
  startsLinear *= stepSpectrum;

  // stop Mel frequencies of filters
  MatrixXR stops(_bandCount, 1);
  for (int i=0; i<stops.rows(); i++) {
    stops(i, 0) = (Real(i + 2) * stepMel + lowMel);
  }

  MatrixXR stopsLinear;
  _melToLinearMatrix(stops, &stopsLinear);
  stopsLinear *= stepSpectrum;


  // center Mel frequencies of filters
  MatrixXR centers(_bandCount, 1);
  for (int i=0; i<centers.rows(); i++) {
    centers(i, 0) = (Real(i + 1) * stepMel + lowMel);
  }

  _centersLinear = startsLinear + (stopsLinear - startsLinear) / 2.0;
  //melToLinearMatrixGreenwood1990(centers, &centersLinear);
  //centersLinear *= stepSpectrum;
  
  // start bins of filters
  MatrixXI startBins = startsLinear.cwise().ceil().cast<Integer>();

  // stop bins of filters
  MatrixXI stopBins = stopsLinear.cwise().ceil().cast<Integer>();

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
      Real center = _centersLinear(i, 0);
      Real stop = stopsLinear(i, 0);
      
      triangleWindow(&newFilter, start - startBin, stop  - startBin, center  - startBin);
    }
    
    weights.push_back(newFilter);
  }

  _bands.setStartsWeights(startBins, weights);

  DEBUG("MELBANDS: Finished set up...");
}

void MelBands::process(const MatrixXR& spectrum, MatrixXR* bands) {  
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

void MelBands::reset(){
  // Initial values
  _bands.reset();
}

void MelBands::starts(MatrixXI* result) const {
  return _bands.starts( result );
}

void MelBands::centers(MatrixXR* result) const {
  (*result) = _centersLinear;
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

Real MelBands::lowFrequency() const{
  return _lowFrequency;
}
  
void MelBands::setLowFrequency( Real frequency, bool callSetup ){
  _lowFrequency = frequency;
  if ( callSetup ) setup();
}

Real MelBands::highFrequency() const{
  return _highFrequency;
}
  
void MelBands::setHighFrequency( Real frequency, bool callSetup ){
  _highFrequency = frequency;
  if ( callSetup ) setup();
}

Real MelBands::samplerate() const{
  return _samplerate;
}
  
void MelBands::setSamplerate( Real frequency, bool callSetup ){
  _samplerate = frequency;
  if ( callSetup ) setup();
}

int MelBands::bandCount() const {
  return _bandCount;
}

void MelBands::setBandCount( int count, bool callSetup ) {
  _bandCount = count;
  if ( callSetup ) setup();
}

int MelBands::fftSize() const{
  return _fftSize;
}

void MelBands::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

MelBands::ScaleType MelBands::scaleType() const{
  return _scaleType;
}

void MelBands::setScaleType( MelBands::ScaleType type, bool callSetup ) {
  _scaleType = type;
  if ( callSetup ) setup();
}
