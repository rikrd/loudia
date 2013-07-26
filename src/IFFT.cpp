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

#include "IFFT.h"

using namespace std;
using namespace Eigen;

IFFT::IFFT(int fftSize, bool zeroPhase) :
  _in( NULL ),
  _out( NULL ),
  _fftplan( NULL )
{
  LOUDIA_DEBUG("IFFT: Constructor fftSize: " << fftSize
        << ", zeroPhase: " << zeroPhase);

  setFftSize( fftSize, false );
  setZeroPhase( zeroPhase, false );
  
  setup();
  
  LOUDIA_DEBUG("IFFT: Constructed");
}

IFFT::~IFFT(){
  if ( _fftplan ) {
    LOUDIA_DEBUG("FFT: Destroying plan");
    fftwf_destroy_plan( _fftplan );
  }

  if ( _in ) {
    LOUDIA_DEBUG("FFT: Destroying in");
    fftwf_free( _in ); 
  }

  if ( _out ) {
    LOUDIA_DEBUG("FFT: Destroying out");
    fftwf_free( _out );
  }
}

void IFFT::setup(){
  LOUDIA_DEBUG("IFFT: Setting up...");

  // Free the ressources if needed 
  // before setting them up
  if ( _fftplan ) {
    LOUDIA_DEBUG("FFT: Destroying plan");
    fftwf_destroy_plan( _fftplan );
  }

  if ( _in ) {
    LOUDIA_DEBUG("FFT: Destroying in");
    fftwf_free( _in ); 
  }

  if ( _out ) {
    LOUDIA_DEBUG("FFT: Destroying out");
    fftwf_free( _out );
  }
  
  _halfSize = ( _fftSize / 2 ) + 1;

  // Allocate the ressources needed  
  _in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _halfSize );
  _out = (Real*) fftwf_malloc(sizeof(Real) * _fftSize);
  
  _fftplan = fftwf_plan_dft_c2r_1d( _fftSize, _in, _out,
                                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT );
      
  LOUDIA_DEBUG("IFFT: Finished set up...");
}

void IFFT::process(const MatrixXC& ffts, MatrixXR* frames){
  const int rows = ffts.rows();
  
  (*frames).resize(rows, _fftSize);

  for (int i = 0; i < rows; i++){    
    // Fill the buffer with zeros
    Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _halfSize) = ffts.row(i);
    
    // Process the data
    fftwf_execute(_fftplan);

    // Take the data from _out
    if(_zeroPhase){

      int half_plus = (int)ceil((Real)_fftSize / 2.0);
      int half_minus = (int)floor((Real)_fftSize / 2.0);
      
      // Take second half of the frame from the beginning 
      (*frames).row(i).block(0, half_minus, 1, half_plus) = Eigen::Map<MatrixXR>(_out, 1, _fftSize).block(0, 0, 1, half_plus) / _fftSize;
      
      // and first half of the frame from the end
      (*frames).row(i).block(0, 0, 1, half_minus) = Eigen::Map<MatrixXR>(_out, 1, _fftSize).block(0, _fftSize - half_minus, 1, half_minus) / _fftSize;

    }else{
      
      // Take all of the frame from the beginning
      (*frames).row(i) = Eigen::Map<MatrixXR>(_out, 1, _fftSize) / _fftSize;

    }
  }
}

void IFFT::reset(){
}

int IFFT::fftSize() const{
  return _fftSize;
}

void IFFT::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

bool IFFT::zeroPhase() const{
  return _zeroPhase;
}

void IFFT::setZeroPhase( bool zeroPhase, bool callSetup ) {
  _zeroPhase = zeroPhase;
  if ( callSetup ) setup();
}
