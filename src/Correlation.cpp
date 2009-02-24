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

Correlation::Correlation(int inputLengthA, int inputLengthB, int maxLag, int minLag) :
  _inputLengthA( inputLengthA ),
  _inputLengthB( inputLengthB ),
  _minLag( max(-max(_inputLengthA, _inputLengthB) + 1, minLag) ),
  _maxLag( min( min(_inputLengthA, _inputLengthB), maxLag) ),
  _useFFT( (_maxLag - _minLag) > 128 ),
  _fftSize(nextPowerOf2(((_inputLengthA + _inputLengthB) - 1)*2)),
  _fft( _fftSize, false ),
  _ifft( _fftSize, false )

{
  DEBUG("CORRELATION: Construction inputLengthA: " << _inputLengthA
        << " inputLengthB: " << _inputLengthB
        << " minLag: " << _minLag
        << " maxLag: " << _maxLag
        << " useFFT: " << _useFFT
        << " fftLength: " << _fftSize);

  setup();
}


Correlation::Correlation(int inputLengthA, int inputLengthB, int maxLag, int minLag, bool useFFT) :
  _inputLengthA( inputLengthA ),
  _inputLengthB( inputLengthB ),
  _minLag( max(-max(_inputLengthA, _inputLengthB) + 1, minLag) ),
  _maxLag( min( min(_inputLengthA, _inputLengthB), maxLag) ),
  _useFFT( useFFT ),
  _fftSize( nextPowerOf2(((_inputLengthA + _inputLengthB) - 1)*2) ),
  _fft( _fftSize, false ),
  _ifft( _fftSize, false )
{
  DEBUG("CORRELATION: Construction inputLengthA: " << _inputLengthA
        << " inputLengthB: " << _inputLengthB
        << " minLag: " << _minLag
        << " maxLag: " << _maxLag
        << " useFFT: " << _useFFT
        << " fftLength: " << _fftSize);
  
  setup();
}

Correlation::~Correlation(){}

void Correlation::setup(){
  // Prepare the buffers
  DEBUG("CORRELATION: Setting up...");

  if ( _useFFT ) {
    _fft.setup();
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
