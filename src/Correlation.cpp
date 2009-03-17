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

#include "Correlation.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

Correlation::Correlation(int inputSizeA, int inputSizeB, int maxLag, int minLag) :
  _inputSizeA( inputSizeA ),
  _inputSizeB( inputSizeB ),
  _minLag( max(-max(_inputSizeA, _inputSizeB) + 1, minLag) ),
  _maxLag( min( min(_inputSizeA, _inputSizeB), maxLag) ),
  _useFFT( (_maxLag - _minLag) > 128 ),
  _fftSize(nextPowerOf2(((_inputSizeA + _inputSizeB) - 1)*2)),
  _fft( _fftSize, false ),
  _ifft( _fftSize, false )

{
  DEBUG("CORRELATION: Construction inputSizeA: " << _inputSizeA
        << " inputSizeB: " << _inputSizeB
        << " minLag: " << _minLag
        << " maxLag: " << _maxLag
        << " useFFT: " << _useFFT
        << " fftSize: " << _fftSize);

  setup();
}


Correlation::Correlation(int inputSizeA, int inputSizeB, int maxLag, int minLag, bool useFFT) :
  _inputSizeA( inputSizeA ),
  _inputSizeB( inputSizeB ),
  _minLag( max(-max(_inputSizeA, _inputSizeB) + 1, minLag) ),
  _maxLag( min( min(_inputSizeA, _inputSizeB), maxLag) ),
  _useFFT( useFFT ),
  _fftSize( nextPowerOf2(((_inputSizeA + _inputSizeB) - 1)*2) ),
  _fft( _fftSize, false ),
  _ifft( _fftSize, false )
{
  DEBUG("CORRELATION: Construction inputSizeA: " << _inputSizeA
        << " inputSizeB: " << _inputSizeB
        << " minLag: " << _minLag
        << " maxLag: " << _maxLag
        << " useFFT: " << _useFFT
        << " fftSize: " << _fftSize);
  
  setup();
}

Correlation::~Correlation(){}

void Correlation::setup(){
  // Prepare the buffers
  DEBUG("CORRELATION: Setting up...");

  if ( _useFFT ) {
    _fftSize = nextPowerOf2( ( (_inputSizeA + _inputSizeB) - 1) * 2 );

    _fft.setFftSize( _fftSize );
    _fft.setup();
    _ifft.setFftSize( _fftSize );
    _ifft.setup();
  }
  
  reset();
  
  DEBUG("CORRELATION: Finished setup.");
}

void Correlation::process(const MatrixXR& inputA, const MatrixXR& inputB, MatrixXR* correlation){
  const int rows = inputA.rows();

  if ( rows != inputB.rows() ) {
    // Thorw ValueError rows of A and B must be the same
  }
  
  (*correlation).resize(rows, _maxLag - _minLag);
  
  if ( _useFFT ) {
    
    _fft.process(inputA, &_fftA);
    _fft.process(inputB, &_fftB);
        
    _ifft.process(_fftA.cwise() * _fftB.conjugate(), &_result);
    
    // TODO: use Eigen rowwise().shift(_fftSize - 2) when it will exist
    for(int i = _minLag;  i < _maxLag; i++ ){
      (*correlation).col(i - _minLag) = _result.col(((_fftSize-2) + (i - _minLag)) % _fftSize);
    }

  } else {
    
    correlate(inputA, inputB, correlation, _minLag, _maxLag);
    
  }
}

void Correlation::reset(){
  // Initial values
}

int Correlation::inputSizeA() const {
  return _inputSizeA;
}
  
void Correlation::setInputSizeA( int size ) {
  _inputSizeA = size;
  setup();
}

int Correlation::inputSizeB() const {
  return _inputSizeB;
}
  
void Correlation::setInputSizeB( int size ) {
  _inputSizeB = size;
  setup();
}

int Correlation::minLag() const {
  return _minLag;
}
  
void Correlation::setMinLag( int lag ) {
  if ( lag >= _maxLag ) {
    // Thorw ValueError, "The minLag should be smaller than the maxLag."
  }

  _minLag = max(-max(_inputSizeA, _inputSizeB) + 1, lag);
  _minLag = min( min(_inputSizeA, _inputSizeB), lag);  
}

int Correlation::maxLag() const {
  return _maxLag;
}
  
void Correlation::setMaxLag( int lag ) {
  if ( lag <= _minLag ) {
    // Thorw ValueError, "The maxLag should be larger than the minLag."
  }
 
  _maxLag = max(-max(_inputSizeA, _inputSizeB) + 1, lag);
  _maxLag = min( min(_inputSizeA, _inputSizeB), lag);  
}

bool Correlation::useFFT() const {
  return _useFFT;
}  

void Correlation::setUseFFT( bool useFFT ) {
  _useFFT = useFFT;
  setup();
}
