/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com.com>
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

#include <cmath>

#include "odfcomplex.h"
#include "window.h"
#include "fft.h"
#include "unwrap.h"

using namespace std;

// import most common Eigen types 
using namespace Eigen;

ODFComplex::ODFComplex(int frameLength, 
                       int fftLength, 
                       Window::WindowType windowType, 
                       bool zeroPhase) : _window(frameLength, windowType), 
                                         _fft(frameLength, fftLength, zeroPhase),
                                         _unwrap((int)(fftLength / 2.0)) {
  
  DEBUG("ODFComplex: Constructor frameLength: " << frameLength << ", fftLength: " << fftLength);
  
  _frameLength = frameLength;
  _fftLength = fftLength;
  _windowType = windowType;
  _zeroPhase = zeroPhase;
  
  setup();
}

ODFComplex::~ODFComplex() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void ODFComplex::setup(){
  // Prepare the buffers
  DEBUG("ODFComplex: Setting up...");

  _window.setup();
  _fft.setup();
  _unwrap.setup();

  reset();
  DEBUG("ODFComplex: Finished set up...");
}


void ODFComplex::process(MatrixXR samples, MatrixXR* odfValue){
  DEBUG("ODFComplex: Processing windowed");

  (*odfValue).resize(1, 1);
  _spectrum.resize(samples.rows(), (int)ceil(_fftLength / 2.0));

  _window.process(samples, &_windowed);

  _fft.process(_windowed, &_ffted);
  
  _spectrum = _ffted.block(0, 0, _ffted.rows(), (int)ceil(_fftLength / 2.0));

  _unwrap.process(_spectrum.angle().real().cast<Real>(), &_unwrappedAngle);

  DEBUG("ODFComplex: Processing unwrappedAngle");
  cout << _unwrappedAngle << endl;
  
  //(*odfValue)(0, 0) = 0.0;
  DEBUG("ODFComplex: Finished Processing");
}

void ODFComplex::reset(){
  // Initial values
  _window.reset();
  _fft.reset();
}
