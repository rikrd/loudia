/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
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

#include "FFT.h"

using namespace std;
using namespace Eigen;

FFT::FFT(int fftSize, bool zeroPhase) :
  _in( NULL ),
  _out( NULL ),
  _fftplan( NULL )
{
  DEBUG("FFT: Constructor fftSize: " << fftSize 
        << ", zeroPhase: " << zeroPhase);

  setFftSize( fftSize, false );
  setZeroPhase( zeroPhase, false );
  
  setup();
  
  DEBUG("FFT: Constructed");
}

FFT::~FFT(){
  DEBUG("FFT: Destroying...");
  if ( _fftplan ) {
    DEBUG("FFT: Destroying plan");
    fftwf_destroy_plan( _fftplan );
  }

  if ( _in ) {
    DEBUG("FFT: Destroying in");
    fftwf_free( _in ); 
  }

  if ( _out ) {
    DEBUG("FFT: Destroying out");
    fftwf_free( _out );
  }
  DEBUG("FFT: Destroyed out");
}

void FFT::setup(){
  DEBUG("FFT: Setting up...");
  
  // Free the ressources if needed 
  // before setting them up
  if ( _fftplan ) {
    DEBUG("FFT: Destroying plan");
    fftwf_destroy_plan( _fftplan );
  }

  if ( _in ) {
    DEBUG("FFT: Destroying in");
    fftwf_free( _in ); 
  }

  if ( _out ) {
    DEBUG("FFT: Destroying out");
    fftwf_free( _out );
  }
  
  _halfSize = ( _fftSize / 2 ) + 1;
  
  // Allocate the ressources needed
  _in = (Real*) fftwf_malloc(sizeof(Real) * _fftSize);
  _out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _halfSize);
  
  _fftplan = fftwf_plan_dft_r2c_1d( _fftSize, _in, _out,
                                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT );
      
  DEBUG("FFT: Finished set up...");
}

void FFT::process(const MatrixXR& frames, MatrixXC* ffts){
  const int cols = frames.cols();
  const int rows = frames.rows();

  if(_fftSize < cols){
    // Throw exception, the FFT size must be greater or equal than the input size
  }
  
  (*ffts).resize(rows, _halfSize);

  for (int i = 0; i < rows; i++){    
    // Fill the buffer with zeros
    Eigen::Map<MatrixXR>(_in, 1, _fftSize) = MatrixXR::Zero(1, _fftSize);
    
    // Put the data in _in
    if(_zeroPhase){

      int half_plus = (int)ceil((Real)cols / 2.0);
      int half_minus = (int)floor((Real)cols / 2.0);
      
      // Put second half of the frame at the beginning 
      Eigen::Map<MatrixXR>(_in, 1, _fftSize).block(0, 0, 1, half_plus) = frames.row(i).block(0, half_minus, 1, half_plus);
      
      // and first half of the frame at the end
      Eigen::Map<MatrixXR>(_in, 1, _fftSize).block(0, _fftSize - half_minus, 1, half_minus) = frames.row(i).block(0, 0, 1, half_minus);


    }else{
      // Put all of the frame at the beginning
      Eigen::Map<MatrixXR>(_in, 1, _fftSize).block(0, 0, 1, cols) = frames.row(i);
    }

    // Process the data
    fftwf_execute(_fftplan);

    // Take the data from _out
    (*ffts).row(i) = Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_out), 1, _halfSize);
  }
}

void FFT::reset(){
}

int FFT::fftSize() const{
  return _fftSize;
}

void FFT::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

bool FFT::zeroPhase() const{
  return _zeroPhase;
}

void FFT::setZeroPhase( bool zeroPhase, bool callSetup ) {
  _zeroPhase = zeroPhase;
  if ( callSetup ) setup();
}
