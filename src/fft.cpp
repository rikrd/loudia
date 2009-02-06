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

#include "typedefs.h"
#include "debug.h"

#include "fft.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

FFT::FFT(int frameSize, int fftSize, bool zeroPhase) :
  _frameSize( frameSize ),
  _fftSize( fftSize ),
  _zeroPhase( zeroPhase ),
  _halfSize( fftSize / 2 + 1 )
{
  DEBUG("FFT: Constructor frameSize: " << frameSize 
        << ", fftSize: " << fftSize 
        << ", zeroPhase: " << zeroPhase);

  if(_fftSize < _frameSize){
    // Throw exception, the FFT size must be greater or equal than the input size
  }
  
  setup();
  
  DEBUG("FFT: Constructed");
}

FFT::~FFT(){
  DEBUG("FFT: Destroying...");
  fftwf_destroy_plan( _fftplan );
  DEBUG("FFT: Destroyed plan");
  fftwf_free( _in ); 
  DEBUG("FFT: Destroyed in");

  DEBUG("FFT: Destroying out");
  fftwf_free( _out );
  DEBUG("FFT: Destroyed out");
}

void FFT::setup(){
  DEBUG("FFT: Setting up...");
  
  _in = (Real*) fftwf_malloc(sizeof(Real) * _fftSize);
  _out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _halfSize);
  
  _fftplan = fftwf_plan_dft_r2c_1d( _fftSize, _in, _out,
                                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT );
      
  DEBUG("FFT: Finished set up...");
}

void FFT::process(const MatrixXR& frames, MatrixXC* ffts){
  (*ffts).resize(frames.rows(), _halfSize);

  for (int i = 0; i < frames.rows(); i++){    
    // Fill the buffer with zeros
    Eigen::Map<MatrixXR>(_in, 1, _fftSize) = MatrixXR::Zero(1, _fftSize);
    
    // Put the data in _in
    if(_zeroPhase){

      int half_plus = ceil((Real)_frameSize / 2.0);
      int half_minus = floor((Real)_frameSize / 2.0);
      
      // Put second half of the frame at the beginning 
      Eigen::Map<MatrixXR>(_in, 1, _fftSize).block(0, 0, 1, half_plus) = frames.row(i).block(0, half_minus, 1, half_plus);
      
      // and first half of the frame at the end
      Eigen::Map<MatrixXR>(_in, 1, _fftSize).block(0, _fftSize - half_minus, 1, half_minus) = frames.row(i).block(0, 0, 1, half_minus);


    }else{

      // Put all of the frame at the beginning
      Eigen::Map<MatrixXR>(_in, 1, _fftSize).block(0, 0, 1, _frameSize) = frames.row(i);
    }
    // Process the data
    fftwf_execute(_fftplan);

    // Take the data from _out
    (*ffts).row(i) = Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_out), 1, _halfSize);
  }
}

void FFT::reset(){
}

int FFT::frameSize() const{
  return _frameSize;
}

int FFT::fftSize() const{
  return _fftSize;
}
