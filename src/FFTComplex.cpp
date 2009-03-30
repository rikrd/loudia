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

#include "FFTComplex.h"

using namespace std;
using namespace Eigen;

FFTComplex::FFTComplex(int frameSize, int fftSize, bool zeroPhase) :
  _in( NULL ),
  _out( NULL ),
  _fftplan( NULL )
{
  DEBUG("FFTComplex: Constructor frameSize: " << frameSize 
        << ", fftSize: " << fftSize 
        << ", zeroPhase: " << zeroPhase);

  if(_fftSize < _frameSize){
    // Throw exception, the FFTComplex size must be greater or equal than the input size
  }

  setFrameSize( frameSize, false );
  setFftSize( fftSize, false );
  setZeroPhase( zeroPhase, false );
  
  setup();
  
  DEBUG("FFTComplex: Constructed");
}

FFTComplex::~FFTComplex(){
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

void FFTComplex::setup(){
  DEBUG("FFTComplex: Setting up...");

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
  
  _in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _fftSize);
  _out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _fftSize);
    
  _fftplan = fftwf_plan_dft_1d( _fftSize, _in, _out,
                                FFTW_FORWARD, FFTW_ESTIMATE );
  
  DEBUG("FFTComplex: Finished set up...");
}

template<typename FrameMatrixType>
void FFTComplex::process(const FrameMatrixType& frames, MatrixXC* ffts){
  (*ffts).resize(frames.rows(), _fftSize);

  for (int i = 0; i < frames.rows(); i++){    
    // Fill the buffer with zeros
    Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _fftSize) = MatrixXC::Zero(1, _fftSize);
    
    // Put the data in _in
    if(_zeroPhase){

      int half_plus = (int)ceil((Real)_frameSize / 2.0);
      int half_minus = (int)floor((Real)_frameSize / 2.0);

      // Put second half of the frame at the beginning 
      Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _fftSize).block(0, 0, 1, half_plus) = frames.row(i).block(0, half_minus, 1, half_plus).template cast<Complex>();
      
      // and first half of the frame at the end
      Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _fftSize).block(0, _fftSize - half_minus, 1, half_minus) = frames.row(i).block(0, 0, 1, half_minus).template cast<Complex>();


    }else{

      // Put all of the frame at the beginning
      Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_in), 1, _fftSize).block(0, 0, 1, _frameSize) = frames.row(i).template cast<Complex>();
    }
    // Process the data
    fftwf_execute(_fftplan);

    // Take the data from _out
    (*ffts).row(i) = Eigen::Map<MatrixXC>(reinterpret_cast< Complex* >(_out), 1, _fftSize);
  }
}

void FFTComplex::process(const MatrixXR& frames, MatrixXC* ffts){
  process<MatrixXR>(frames, ffts);
}

void FFTComplex::process(const MatrixXC& frames, MatrixXC* ffts){
  process<MatrixXC>(frames, ffts);
}

void FFTComplex::reset(){
}

int FFTComplex::fftSize() const{
  return _fftSize;
}

void FFTComplex::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

int FFTComplex::frameSize() const{
  return _frameSize;
}

void FFTComplex::setFrameSize( int size, bool callSetup ) {
  _frameSize = size;
  if ( callSetup ) setup();
}

bool FFTComplex::zeroPhase() const{
  return _zeroPhase;
}

void FFTComplex::setZeroPhase( bool zeroPhase, bool callSetup ) {
  _zeroPhase = zeroPhase;
  if ( callSetup ) setup();
}
