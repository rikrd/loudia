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

Autocorrelation::Autocorrelation(int inputSize, int maxLag, int minLag) :
  _inputSize( inputSize ),
  _minLag( max(-inputSize + 1, minLag) ),
  _maxLag( min(inputSize, maxLag) ),
  _useFFT( (_maxLag - _minLag) > 128 ),
  _fft( nextPowerOf2(((_maxLag - _minLag)-1)*2), false ),
  _ifft( nextPowerOf2(((_maxLag - _minLag)-1)*2), false )
{
  DEBUG("AUTOCORRELATION: Construction inputSize: " << _inputSize
        << " minLag: " << _minLag
        << " maxLag: " << _maxLag
        << " useFFT: " << _useFFT);

  setup();
}


Autocorrelation::Autocorrelation(int inputSize, int maxLag, int minLag, bool useFFT) :
  _inputSize( inputSize ),
  _minLag( max(-inputSize + 1, minLag) ),
  _maxLag( min(inputSize, maxLag) ),
  _useFFT( useFFT ),
  _fft( nextPowerOf2(((_maxLag - _minLag)-1)*2), false ),
  _ifft( nextPowerOf2(((_maxLag - _minLag)-1)*2), false )
{
  DEBUG("AUTOCORRELATION: Construction inputSize: " << _inputSize
        << " minLag: " << _minLag
        << " maxLag: " << _maxLag
        << " useFFT: " << _useFFT);
  
  setup();
}

Autocorrelation::~Autocorrelation(){}

void Autocorrelation::setup(){
  // Prepare the buffers
  DEBUG("AUTOCORRELATION: Setting up...");

  if ( _useFFT ) {
    _fft.setup();
    _ifft.setup();
  }
  
  reset();
  
  DEBUG("AUTOCORRELATION: Finished setup.");
}

void Autocorrelation::process(const MatrixXR& frames, MatrixXR* autocorrelation){
  const int rows = frames.rows();

  (*autocorrelation).resize(rows, _maxLag - _minLag);
  
  if ( _useFFT ) {

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
  
void Autocorrelation::setInputSize( int size ) {
  _inputSize = size;
}

int Autocorrelation::minLag() const {
  return _minLag;
}
  
void Autocorrelation::setMinLag( int lag ) {
  _minLag = max(-_inputSize + 1, lag);  
}

int Autocorrelation::maxLag() const {
  return _maxLag;
}
  
void Autocorrelation::setMaxLag( int lag ) {
  _maxLag = min(_inputSize, lag);
}

bool Autocorrelation::useFFT() const {
  return _useFFT;
}  

void Autocorrelation::setUseFFT( bool useFFT ) {
  _useFFT = useFFT;
}
