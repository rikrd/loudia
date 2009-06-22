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

#include "Correlation.h"

#include "Utils.h"

using namespace std;
using namespace Eigen;

Correlation::Correlation()
{
  LOUDIA_DEBUG("CORRELATION: Construction inputSizeA: " << inputSizeA
        << " inputSizeB: " << inputSizeB
        << " minLag: " << minLag
        << " maxLag: " << maxLag);

  setInputSizeA( 1024, false );
  setInputSizeB( 1024, false );
  setMinLag( -std::numeric_limits<int>::max(), false );
  setMaxLag( std::numeric_limits<int>::max(), false );
  setUseFft( (_maxLag - _minLag) > 128, false );

  setup();
}

Correlation::Correlation(int inputSizeA, int inputSizeB, int maxLag, int minLag)
{
  LOUDIA_DEBUG("CORRELATION: Construction inputSizeA: " << inputSizeA
        << " inputSizeB: " << inputSizeB
        << " minLag: " << minLag
        << " maxLag: " << maxLag);

  setInputSizeA( inputSizeA, false );
  setInputSizeB( inputSizeB, false );
  setMinLag( minLag, false );
  setMaxLag( maxLag, false );
  setUseFft( (_maxLag - _minLag) > 128, false );

  setup();
}


Correlation::Correlation(int inputSizeA, int inputSizeB, int maxLag, int minLag, bool useFft)
{
  LOUDIA_DEBUG("CORRELATION: Construction inputSizeA: " << inputSizeA
        << " inputSizeB: " << inputSizeB
        << " minLag: " << minLag
        << " maxLag: " << maxLag
        << " useFft: " << useFft);
  
  setInputSizeA( inputSizeA, false );
  setInputSizeB( inputSizeB, false );
  setMinLag( minLag, false );
  setMaxLag( maxLag, false );
  setUseFft( useFft, false );
  
  setup();
}

Correlation::~Correlation(){}

void Correlation::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("CORRELATION: Setting up...");

  if ( _useFft ) {
    _fftSize = nextPowerOf2( ( (_inputSizeA + _inputSizeB) - 1) * 2 );
    
    _fft.setFftSize( _fftSize, false );
    _fft.setZeroPhase( false, false );
    _fft.setup();

    _ifft.setFftSize( _fftSize, false );
    _ifft.setZeroPhase( false, false );
    _ifft.setup();
  }
  
  reset();
  
  LOUDIA_DEBUG("CORRELATION: Finished setup.");
}

void Correlation::process(const MatrixXR& inputA, const MatrixXR& inputB, MatrixXR* correlation){
  const int rows = inputA.rows();

  if ( rows != inputB.rows() ) {
    // Thorw ValueError rows of A and B must be the same
  }
  
  (*correlation).resize(rows, _maxLag - _minLag);
  
  if ( _useFft ) {
    
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
  
void Correlation::setInputSizeA( int size, bool callSetup ) {
  _inputSizeA = size;
  if ( callSetup ) setup();
}

int Correlation::inputSizeB() const {
  return _inputSizeB;
}
  
void Correlation::setInputSizeB( int size, bool callSetup ) {
  _inputSizeB = size;
  if ( callSetup ) setup();
}

int Correlation::minLag() const {
  return _minLag;
}
  
void Correlation::setMinLag( int lag, bool callSetup ) {
  if ( lag >= _maxLag ) {
    // Thorw ValueError, "The minLag should be smaller than the maxLag."
  }

  _minLag = max(-max(_inputSizeA, _inputSizeB) + 1, lag);
  _minLag = min( min(_inputSizeA, _inputSizeB), _minLag);  
  if ( callSetup ) setup();
}

int Correlation::maxLag() const {
  return _maxLag;
}
  
void Correlation::setMaxLag( int lag, bool callSetup ) {
  if ( lag <= _minLag ) {
    // Thorw ValueError, "The maxLag should be larger than the minLag."
  }
 
  _maxLag = max(-max(_inputSizeA, _inputSizeB) + 1, lag);
  _maxLag = min( min(_inputSizeA, _inputSizeB), _maxLag);  
  if ( callSetup ) setup();
}

bool Correlation::useFft() const {
  return _useFft;
}  

void Correlation::setUseFft( bool useFft, bool callSetup ) {
  _useFft = useFft;
  if ( callSetup ) setup();
}
