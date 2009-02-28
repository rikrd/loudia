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

SpectralWhitening::SpectralWhitening(int fftSize, Real f0, Real f1, Real samplerate, Real compressionFactor, int numBands, MelBands::ScaleType scaleType) :
  _fftSize( fftSize ),
  _halfSize( ( _fftSize / 2 ) + 1 ),
  _f0( f0 ),
  _f1( f1 ),
  _samplerate( samplerate ),
  _compressionFactor( compressionFactor ),
  _numBands( numBands ),
  _scaleType( scaleType ),
  _bands( _f0, _f1, _numBands, _samplerate, _fftSize, _scaleType)
{
  DEBUG("SPECTRALWHITENING: Construction fftSize: " << _fftSize
        << " samplerate: " << _samplerate
        << " compressionFactor: " << _compressionFactor
        << " numBands: " << _numBands
        << " f0: " << _f0
        << " f1: " << _f1 );

  setup();
}

SpectralWhitening::~SpectralWhitening(){}

void SpectralWhitening::setup(){
  DEBUG("SPECTRALWHITENING: Setting up...");

  // Setup the bands
  _bands.setup();

  _bands.centers(&_centers);
  
  cout << _centers << endl;

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
  for ( int band = 1; band < _numBands; band++ ) {
    for (; col < _centers(band, 0); col++ ) {
      _compressionWeights.col(col) = (((Real)col - _centers(band - 1, 0)) * (_bandEnergy.col(band) - _bandEnergy.col(band-1)) / (_centers(band, 0) - _centers(band - 1, 0))) + _bandEnergy.col(band - 1);
    }
  }

  // Interpolate the region after the last frequency center
  for (; col < _halfSize; col++ ) {
      _compressionWeights.col(col) = (((Real)col - _centers(_numBands - 1, 0)) * ( -_bandEnergy.col(_numBands - 1)) / (_halfSize - _centers(_numBands - 1, 0))) + _bandEnergy.col(_numBands - 1);
  }

  // Apply compression weihgts
  (*result) = spectrum.cwise() * _compressionWeights;
}

void SpectralWhitening::reset(){
  // Initial values

  _bands.reset();
}
