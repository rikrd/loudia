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
  _bands( _f0, _f1, _numBands, _samplerate, _fftSize, _scaleType),
  _resample( _numBands, _halfSize, (Real)_halfSize / (Real)_numBands )
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

  reset();

  DEBUG("SPECTRALWHITENING: Finished setup.");
}

void SpectralWhitening::process(const MatrixXR& spectrum, MatrixXR* result){
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();

  (*result).resize( rows, cols );

  _bands.process( spectrum.cwise().square(), &_bandEnergy );
  
  _resample.process( (_bandEnergy / _fftSize).cwise().sqrt().cwise().pow(_compressionFactor - 1.0), &_compressionWeights );

  (*result) = spectrum.cwise() * _compressionWeights;
}

void SpectralWhitening::reset(){
  // Initial values

  _bands.reset();
}
