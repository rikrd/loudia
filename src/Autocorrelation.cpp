/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 3 of the License, or   
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

Autocorrelation::Autocorrelation(int inputSize)
{
  LOUDIA_DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize);

  setInputSize( inputSize, false );
  setMinLag( 0, false );
  setMaxLag( inputSize, false );
  setUseFft( (_maxLag - _minLag) > 128, false );

  setup();
}

Autocorrelation::Autocorrelation(int inputSize, int maxLag)
{
  LOUDIA_DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize
        << " maxLag: " << maxLag);

  setInputSize( inputSize, false );
  setMinLag( 0, false );
  setMaxLag( maxLag, false );
  setUseFft( (_maxLag - _minLag) > 128, false );

  setup();
}


Autocorrelation::Autocorrelation(int inputSize, int maxLag, int minLag)
{
  LOUDIA_DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize
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
  LOUDIA_DEBUG("AUTOCORRELATION: Construction inputSize: " << inputSize
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
  LOUDIA_DEBUG("AUTOCORRELATION: Setting up...");

  _calcMinLag = min(_inputSize + 1, max(-_inputSize + 1, _minLag));  
  _calcMaxLag = min(_inputSize + 1, max(-_inputSize + 1, _maxLag));

  if ( _useFft ) {
    _fft.setFftSize( nextPowerOf2( _inputSize * 2, false ) );
    _fft.setZeroPhase( false, false );
    
    _ifft.setFftSize( nextPowerOf2( _inputSize * 2, false ) );
    _ifft.setZeroPhase( false, false );

    _fft.setup();
    _ifft.setup();
  }
  
  reset();
  
  LOUDIA_DEBUG("AUTOCORRELATION: Finished setup.");
}

void Autocorrelation::process(const MatrixXR& frames, MatrixXR* autocorrelation){
  const int rows = frames.rows();

  (*autocorrelation).resize(rows, _maxLag - _minLag);
  
  (*autocorrelation).setZero();
  
  if ( _useFft ) {
    _fft.process(frames, &_tempFft);
    
    _tempFft.cwise() *= _tempFft.conjugate();
    
    _ifft.process(_tempFft, &_temp);

    (*autocorrelation).block(0, _calcMinLag - _minLag, rows, _calcMaxLag - _calcMinLag) = _temp.block(0, 0, rows, _calcMaxLag - _calcMinLag);

  } else {
    correlate(frames, frames, &_temp, _calcMinLag, _calcMaxLag);
    
    (*autocorrelation).block(0, _calcMinLag - _minLag, rows, _calcMaxLag - _calcMinLag) = _temp;
    
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
  _minLag = lag;
  if ( callSetup ) setup();
}

int Autocorrelation::maxLag() const {
  return _maxLag;
}
  
void Autocorrelation::setMaxLag( int lag, bool callSetup ) {
  _maxLag = lag;
  if ( callSetup ) setup();
}

bool Autocorrelation::useFft() const {
  return _useFft;
}  

void Autocorrelation::setUseFft( bool useFft, bool callSetup ) {
  _useFft = useFft;
  if ( callSetup ) setup();
}
