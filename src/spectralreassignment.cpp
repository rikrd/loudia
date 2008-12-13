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

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>
#include <cmath>

#include "window.h"
#include "spectralreassignment.h"
#include "debug.h"
#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

SpectralReassignment::SpectralReassignment(int frameSize, int fftSize, Real samplerate, Window::WindowType windowType) : _windowAlgo(frameSize, windowType), _windowIntegAlgo(frameSize, Window::CUSTOM), _windowDerivAlgo(frameSize, Window::CUSTOM), _fftAlgo(frameSize, fftSize){
  DEBUG("SPECTRALREASSIGNMENT: Constructor frameSize: " << frameSize << ", fftSize: " << fftSize << ", samplerate: " << samplerate << ", windowType: " << windowType);

  _frameSize = frameSize;
  _fftSize = fftSize;
  _samplerate = samplerate;
  _windowType = windowType;
  
  DEBUG("SPECTRALREASSIGNMENT: Constructed");
}

SpectralReassignment::~SpectralReassignment(){}

void SpectralReassignment::setup(){
  DEBUG("SPECTRALREASSIGNMENT: Setting up...");
  
  // Create the time vector
  DEBUG("SPECTRALREASSIGNMENT: Creating time vector...");
  Real timestep = 1.0 / _samplerate;
  _time.resize(_frameSize, _fftSize);
  for(int i = 0; i < _time.rows(); i++){
    _time.row(i) = MatrixXR::Constant(1, _time.cols(), timestep).row(0) * i;
  }

  // Create the freq vector
  DEBUG("SPECTRALREASSIGNMENT: Creating freq vector...");
  Real freqstep = 2.0 * M_PI / _fftSize;
  _freq.resize(_frameSize, _fftSize);
  for(int i = 0; i < _freq.cols(); i++){
    _freq.col(i) = MatrixXR::Constant(_freq.rows(), 1, freqstep).col(0) * i;
  }
  
  // Create the reassign operator matrix
  DEBUG("SPECTRALREASSIGNMENT: Creating reassignment op...");
  _reassignOp.resize(_frameSize, _fftSize);

  // Calculate and set the time weighted window
  DEBUG("SPECTRALREASSIGNMENT: Calculate time weighted window...");
  MatrixXR windowInteg = _windowIntegAlgo.window();
  windowInteg = windowInteg.cwise() * _time.col(0).transpose();
  _windowIntegAlgo.setWindow(windowInteg);

  // Calculate and set the time derivated window
  DEBUG("SPECTRALREASSIGNMENT: Calculate time derivative window...");
  MatrixXR windowDeriv = _windowDerivAlgo.window();
  for(int i = windowDeriv.cols() - 1; i > 0; i--){
    windowDeriv(0, i) = (windowDeriv(0, i) - windowDeriv(0, i-1)) / timestep;
  }
  // TODO: Check what is the initial condition for the window
  // Should this be 0 or just the value it was originally * dt
  windowDeriv(0, 0) = 0.0;
  _windowDerivAlgo.setWindow(windowDeriv);

  // Create the necessary buffers for the windowing
  _window.resize(1, _frameSize);
  _windowInteg.resize(1, _frameSize);
  _windowDeriv.resize(1, _frameSize);

  // Create the necessary buffers for the FFT
  _fft.resize(1, _fftSize);
  _fftInteg.resize(1, _fftSize);
  _fftDeriv.resize(1, _fftSize);
  
  // Setup the algos
  _windowAlgo.setup();
  _windowIntegAlgo.setup();
  _windowDerivAlgo.setup();
  _fftAlgo.setup();

  DEBUG("SPECTRALREASSIGNMENT: Finished set up...");
}


template<class F, class W>
void SpectralReassignment::process(F frames, W* reassigned, W* fft){
  for (int i = 0; i < frames.rows(); i++){

    // Process the windowing
    _windowAlgo.process(frames, &_window);
    _windowIntegAlgo.process(frames, &_windowInteg);
    _windowDerivAlgo.process(frames, &_windowDeriv);

    // Process the FFT
    _fftAlgo.process(_window, &_fft);
    _fftAlgo.process(_windowInteg, &_fftInteg);
    _fftAlgo.process(_windowDeriv, &_fftDeriv);

    // Reassign
    
  }
}

void SpectralReassignment::process(MatrixXC frames, MatrixXC* reassigned, MatrixXC* fft){
  process<MatrixXC, MatrixXC>(frames, reassigned, fft);
}

void SpectralReassignment::process(MatrixXR frames, MatrixXC* reassigned, MatrixXC* fft){
  process<MatrixXR, MatrixXC>(frames, reassigned, fft);
}

void SpectralReassignment::process(MatrixXC frames, MatrixXC* reassigned){
  process<MatrixXC, MatrixXC>(frames, reassigned, &_fft);
}

void SpectralReassignment::process(MatrixXR frames, MatrixXC* reassigned){
  process<MatrixXR, MatrixXC>(frames, reassigned, &_fft);
}

void SpectralReassignment::reset(){
  _windowAlgo.reset();
  _windowIntegAlgo.reset();
  _windowDerivAlgo.reset();
  _fftAlgo.reset();
}

int SpectralReassignment::frameSize() const{
  return _frameSize;
}

int SpectralReassignment::fftSize() const{
  return _fftSize;
}

Window::WindowType SpectralReassignment::windowType() const{
  return _windowType;
}
