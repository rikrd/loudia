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

#include "SpectralNoiseSuppression.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;

SpectralNoiseSuppression::SpectralNoiseSuppression(int fftSize, Real f0, Real f1, Real samplerate) :
  _fftSize( fftSize ),
  _samplerate( samplerate ),
  _k0( f0 / _samplerate ),
  _k1( f1 / _samplerate )
{
  DEBUG("SpectralNoiseSuppression: Construction fftSize: " << _fftSize
        << " samplerate: " << _samplerate
        << " k0: " << _k0
        << " k1: " << _k1);

  setup();
}

SpectralNoiseSuppression::~SpectralNoiseSuppression(){}

void SpectralNoiseSuppression::setup(){
  DEBUG("SPECTRALNOISESUPPRESSION: Setting up...");
  
  // Prepare the bands for the moving average
  int _halfSize = (_fftSize / 2) + 1;
  
  
  _bands.setup()

  reset();

  DEBUG("SPECTRALNOISESUPPRESSION: Finished setup.");
}

void SpectralNoiseSuppression::process(const MatrixXR& spectrum, MatrixXR* result){
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();

  (*result).resize(rows, cols);
  
  (*result) = spectrum;


  // Calculate the wrapped magnitude of the spectrum
  _g = (1 / (_k1 - _k0 + 1) * spectrum.block(0, _k0, rows, _k1 - _k0).cwise().pow(1.0/3.0).rowwise().sum()).cwise().cube();

  for ( int i = 0; i < cols; i++ ) {
    (*result).col(i) = (((*result).col(i).cwise() * _g.inverse()).cwise() + 1.0).cwise().log();
  }

  // Estimate spectral noise
  _bands.process((*result), &_noise);

  // Suppress spectral noise
  (*result) = ((*result) - _noise).cwise().clipUnder();
}

void SpectralNoiseSuppression::reset(){
  // Initial values

  _bands.reset();
}
