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

#include "Autocorrelation.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

Autocorrelation::Autocorrelation(int inputSize, int maxLag, int minLag)
{
  DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize
        << " minLag: " << minLag
        << " maxLag: " << maxLag);

  setInputSize( inputSize, false );
  setMinLag( minLag, false );
  setMaxLag( maxLag, false );
  setUseFft( (_maxLag - _minLag) > 128, false );

  setup();
}


Autocorrelation::Autocorrelation(int inputSize, int maxLag, int minLag, bool useFft)
{
  DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize
        << " minLag: " << minLag
        << " maxLag: " << maxLag
        << " useFft: " << useFft);
  
  setInputSize( inputSize, false );
  setMinLag( minLag, false );
  setMaxLag( maxLag, false );
  setUseFft( useFft, false );
  
  setup();
}

Autocorrelation::~Autocorrelation(){}

void Autocorrelation::setup(){
  // Prepare the buffers
  DEBUG("AUTOCORRELATION: Setting up...");

  if ( _useFft ) {
    _fft.setFftSize( nextPowerOf2(((_maxLag - _minLag)-1)*2), false );
    _fft.setZeroPhase( false, false );
    
    _ifft.setFftSize( nextPowerOf2(((_maxLag - _minLag)-1)*2), false );
    _ifft.setZeroPhase( false, false );

    _fft.setup();
    _ifft.setup();
  }
  
  reset();
  
  DEBUG("AUTOCORRELATION: Finished setup.");
}

void Autocorrelation::process(const MatrixXR& frames, MatrixXR* autocorrelation){
  const int rows = frames.rows();

  (*autocorrelation).resize(rows, _maxLag - _minLag);
  
  if ( _useFft ) {

    MatrixXC temp;
    MatrixXR temp2;
    _fft.process(frames, &temp);
    
    temp.cwise() *= temp.conjugate();
    
    _ifft.process(temp, &temp2);
    
    (*autocorrelation) = temp2.block(0, 0, rows, _maxLag - _minLag);

  } else {
    
    correlate(frames, frames, autocorrelation, _minLag, _maxLag);

  }
}

void Autocorrelation::reset(){
  // Initial values
}

int Autocorrelation::inputSize() const {
  return _inputSize;
}
  
void Autocorrelation::setInputSize( int size, bool callSetup ) {
  _inputSize = size;
  if ( callSetup ) setup();
}

int Autocorrelation::minLag() const {
  return _minLag;
}
  
void Autocorrelation::setMinLag( int lag, bool callSetup ) {
  _minLag = min(_inputSize, max(-_inputSize + 1, lag));  
  if ( callSetup ) setup();
}

int Autocorrelation::maxLag() const {
  return _maxLag;
}
  
void Autocorrelation::setMaxLag( int lag, bool callSetup ) {
  _maxLag = min(_inputSize, max(-_inputSize + 1, lag));
  if ( callSetup ) setup();
}

bool Autocorrelation::useFft() const {
  return _useFft;
}  

void Autocorrelation::setUseFft( bool useFft, bool callSetup ) {
  _useFft = useFft;
  if ( callSetup ) setup();
}
