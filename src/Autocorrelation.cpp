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

Autocorrelation::Autocorrelation(int inputLength, int maxLag, int minLag) :
  _inputLength( inputLength ),
  _minLag( minLag ),
  _maxLag( min(inputLength, maxLag) ),
  _useFFT( (_maxLag - _minLag) > 128 ),
  _fft( inputLength, inputLength, false ),
  _ifft( inputLength, inputLength, false )
{
  DEBUG("AUTOCORRELATION: Construction inputLength: " << _inputLength
        << " minLag: " << _minLag
        << " maxLag: " << _maxLag
        << " useFFT: " << _useFFT);

  setup();
}


Autocorrelation::Autocorrelation(int inputLength, int maxLag, int minLag, bool useFFT) :
  _inputLength( inputLength ),
  _minLag( minLag ),
  _maxLag( min(inputLength, maxLag) ),
  _useFFT( useFFT ),
  _fft( inputLength, inputLength, false ),
  _ifft( inputLength, inputLength, false )
{
  DEBUG("AUTOCORRELATION: Construction inputLength: " << _inputLength
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

void Autocorrelation::process(const MatrixXR& input, MatrixXR* autocorrelation){
  const int rows = input.rows();

  (*autocorrelation).resize(rows, _maxLag - _minLag);
  
  if ( _useFFT ) {

    MatrixXC temp;
    MatrixXR temp2;
    _fft.process(input, &temp);
    
    temp.cwise() *= temp.conjugate();
    
    _ifft.process(temp, &temp2);
    
    (*autocorrelation) = temp2.block(0, 0, rows, _maxLag - _minLag);

  } else {

    correlate(input, input, autocorrelation, _minLag, _maxLag);

  }
}

void Autocorrelation::reset(){
  // Initial values
}
