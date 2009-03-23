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

#include "SpectralWhitening.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

SpectralWhitening::SpectralWhitening(int fftSize, Real lowFrequency, Real highFrequency, Real samplerate, Real compressionFactor, int bandCount, MelBands::ScaleType scaleType) :
  _fftSize( fftSize ),
  _halfSize( ( _fftSize / 2 ) + 1 ),
  _lowFrequency( lowFrequency ),
  _highFrequency( highFrequency ),
  _samplerate( samplerate ),
  _compressionFactor( compressionFactor ),
  _bandCount( bandCount ),
  _scaleType( scaleType ),
  _bands( _lowFrequency, _highFrequency, _bandCount, _samplerate, _fftSize, _scaleType)
{
  DEBUG("SPECTRALWHITENING: Construction fftSize: " << _fftSize
        << " samplerate: " << _samplerate
        << " compressionFactor: " << _compressionFactor
        << " bandCount: " << _bandCount
        << " lowFrequency: " << _lowFrequency
        << " highFrequency: " << _highFrequency );

  setup();
}

SpectralWhitening::~SpectralWhitening(){}

void SpectralWhitening::setup(){
  DEBUG("SPECTRALWHITENING: Setting up...");

  _halfSize = ( _fftSize / 2 ) + 1;
  
  // Setup the bands
  _bands.setLowFrequency( _lowFrequency, false );
  _bands.setHighFrequency( _highFrequency, false );
  _bands.setBandCount(_bandCount, false );
  _bands.setSamplerate( _samplerate, false );
  _bands.setFftSize(_fftSize, false );
  _bands.setScaleType( _scaleType, false ); 
  _bands.setup();

  _bands.centers(&_centers);
  
  reset();

  DEBUG("SPECTRALWHITENING: Finished setup.");
}

void SpectralWhitening::process(const MatrixXR& spectrum, MatrixXR* result){
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();

  (*result).resize( rows, cols );

  _compressionWeights.resize(rows, _halfSize);

  // Calculate the energy per band
  _bands.process( spectrum.cwise().square(), &_bandEnergy );

  // Calculate compress weights of bands
  _bandEnergy = (_bandEnergy / _fftSize).cwise().sqrt().cwise().pow(_compressionFactor - 1.0);

  // Interpolate linearly between the center frequencies of bands
  // Interpolate the region before the first frequency center
  int col = 0;
  for (; col < _centers(0, 0); col++ ) {
    _compressionWeights.col( col ) = (Real)col * _bandEnergy.col(1) / _centers(0, 0);
  }

  // Interpolate the region between the first and last frequency centers
  for ( int band = 1; band < _bandCount; band++ ) {
    for (; col < _centers(band, 0); col++ ) {
      _compressionWeights.col(col) = (((Real)col - _centers(band - 1, 0)) * (_bandEnergy.col(band) - _bandEnergy.col(band-1)) / (_centers(band, 0) - _centers(band - 1, 0))) + _bandEnergy.col(band - 1);
    }
  }

  // Interpolate the region after the last frequency center
  for (; col < _halfSize; col++ ) {
      _compressionWeights.col(col) = (((Real)col - _centers(_bandCount - 1, 0)) * ( -_bandEnergy.col(_bandCount - 1)) / (_halfSize - _centers(_bandCount - 1, 0))) + _bandEnergy.col(_bandCount - 1);
  }

  // Apply compression weihgts
  (*result) = spectrum.cwise() * _compressionWeights;
}

void SpectralWhitening::reset(){
  // Initial values

  _bands.reset();
}

Real SpectralWhitening::lowFrequency() const{
  return _lowFrequency;
}
  
void SpectralWhitening::setLowFrequency( Real frequency, bool callSetup ){
  _lowFrequency = frequency;
  if ( callSetup ) setup();
}

Real SpectralWhitening::highFrequency() const{
  return _highFrequency;
}
  
void SpectralWhitening::setHighFrequency( Real frequency, bool callSetup ){
  _highFrequency = frequency;
  if ( callSetup ) setup();
}

Real SpectralWhitening::samplerate() const{
  return _samplerate;
}
  
void SpectralWhitening::setSamplerate( Real frequency, bool callSetup ){
  _samplerate = frequency;
  if ( callSetup ) setup();
}

int SpectralWhitening::bandCount() const {
  return _bandCount;
}

void SpectralWhitening::setBandCount( int count, bool callSetup ) {
  _bandCount = count;
  if ( callSetup ) setup();
}

int SpectralWhitening::fftSize() const{
  return _fftSize;
}

void SpectralWhitening::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

MelBands::ScaleType SpectralWhitening::scaleType() const{
  return _scaleType;
}

void SpectralWhitening::setScaleType( MelBands::ScaleType type, bool callSetup ) {
  _scaleType = type;
  if ( callSetup ) setup();
}

Real SpectralWhitening::compressionFactor() const{
  return _compressionFactor;
}
  
void SpectralWhitening::setCompressionFactor( Real factor, bool callSetup ){
  _compressionFactor = factor;
  if ( callSetup ) setup();
}
