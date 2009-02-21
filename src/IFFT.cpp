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

#include "IFFT.h"

using namespace std;
using namespace Eigen;

IFFT::IFFT(int fftSize, int frameSize, bool zeroPhase) :
  _fftSize( fftSize ),
  _frameSize( frameSize ),
  _zeroPhase( zeroPhase ),
  _halfSize( fftSize / 2 + 1 )
{
  DEBUG("IFFT: Constructor frameSize: " << frameSize 
        << ", fftSize: " << fftSize 
        << ", zeroPhase: " << zeroPhase);

  if(_fftSize < _frameSize){
    // Throw exception, the IFFT size must be greater or equal than the input size
  }
  
  setup();
  
  DEBUG("IFFT: Constructed");
}

IFFT::~IFFT(){
  fftwf_destroy_plan( _fftplan );
  fftwf_free( _in ); 
  fftwf_free( _out );
}

void IFFT::setup(){
  DEBUG("IFFT: Setting up...");
  
  _in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _halfSize );
  _out = (Real*) fftwf_malloc(sizeof(Real) * _fftSize);
  
  _fftplan = fftwf_plan_dft_c2r_1d( _fftSize, _in, _out,
                                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT );
      
  DEBUG("IFFT: Finished set up...");
}

void IFFT::process(const MatrixXC& ffts, MatrixXR* frames){
  const int rows = ffts.rows();
  
  (*frames).resize(rows, _frameSize);

  for (int i = 0; i < rows; i++){    
    // Fill the buffer with zeros
    Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _halfSize) = ffts.row(i);
    
    // Process the data
    fftwf_execute(_fftplan);

    // Take the data from _out
    if(_zeroPhase){

      int half_plus = ceil((Real)_frameSize / 2.0);
      int half_minus = floor((Real)_frameSize / 2.0);
      
      // Take second half of the frame from the beginning 
      (*frames).row(i).block(0, half_minus, 1, half_plus) = Eigen::Map<MatrixXR>(_out, 1, _fftSize).block(0, 0, 1, half_plus) / _fftSize;
      
      // and first half of the frame from the end
      (*frames).row(i).block(0, 0, 1, half_minus) = Eigen::Map<MatrixXR>(_out, 1, _fftSize).block(0, _fftSize - half_minus, 1, half_minus) / _fftSize;

    }else{

      // Take all of the frame from the beginning
      (*frames).row(i) = Eigen::Map<MatrixXR>(_out, 1, _fftSize).block(0, 0, 1, _frameSize) / _fftSize;
    }
  }
}

void IFFT::reset(){
}

int IFFT::frameSize() const{
  return _frameSize;
}

int IFFT::fftSize() const{
  return _fftSize;
}