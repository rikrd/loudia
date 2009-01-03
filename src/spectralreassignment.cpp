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

#include <Eigen/Core>
#include <Eigen/Array>
#include <iostream>
#include <cmath>

#include "window.h"
#include "spectralreassignment.h"




using namespace std;

// import most common Eigen types 
using namespace Eigen;

SpectralReassignment::SpectralReassignment(int frameSize, int fftSize, Real samplerate, Window::WindowType windowType) : _windowAlgo(frameSize, windowType), _windowIntegAlgo(frameSize, Window::CUSTOM), _windowDerivAlgo(frameSize, Window::CUSTOM), _fftAlgo(frameSize, fftSize, true){
  DEBUG("SPECTRALREASSIGNMENT: Constructor frameSize: " << frameSize << ", fftSize: " << fftSize << ", samplerate: " << samplerate << ", windowType: " << windowType);

  _frameSize = frameSize;
  _fftSize = fftSize;
  _samplerate = samplerate;
  _windowType = windowType;

  setup();
  DEBUG("SPECTRALREASSIGNMENT: Constructed");
}

SpectralReassignment::~SpectralReassignment(){}

void SpectralReassignment::setup(){
  DEBUG("SPECTRALREASSIGNMENT: Setting up...");
  
  // Setup the window so it gets calculated and can be reused
  _windowAlgo.setup();
  
  // Create the time vector
  DEBUG("SPECTRALREASSIGNMENT: Creating time vector...");
  Real timestep = 1.0 / _samplerate;

  // The unit of the vectors is Time Sample fractions
  // So the difference between one coeff and the next is 1
  // and the center of the window must be 0, so even sized windows
  // will have the two center coeffs to -0.5 and 0.5
  // This should be a line going from [-(window_size - 1)/2 ... (window_size - 1)/2]
  _time.resize(_frameSize, 1);
  for(int i = 0; i < _time.rows(); i++){
    _time(i, 0) = (i - Real(_time.rows() - 1)/2.0);
  }

  // Create the freq vector
  DEBUG("SPECTRALREASSIGNMENT: Creating freq vector...");
  Real freqstep = 2.0 * M_PI / _fftSize;

  // The unit of the vectors is Frequency Bin fractions
  // TODO: Must rethink how the frequency vector is initialized
  // as we did for the time vector
  _freq.resize(1, _fftSize);
  for(int i = 0; i < _freq.cols(); i++){
    _freq(0, i) = i;
  }
  
  // Create the reassign operator matrix
  // TODO: check if the current timestamp is enough for a good reassignment
  // we might need for it to depend on past frames (if the reassignment of time
  // goes further than one)
  DEBUG("SPECTRALREASSIGNMENT: Creating reassignment op...");
  _reassignTime.resize(1, _fftSize);
  _reassignFreq.resize(1, _fftSize);

  // Calculate and set the time weighted window
  DEBUG("SPECTRALREASSIGNMENT: Calculate time weighted window...");
  MatrixXR windowInteg = _windowAlgo.window();
  windowInteg = windowInteg.cwise() * _time.transpose();
  _windowIntegAlgo.setWindow(windowInteg);

  // Calculate and set the time derivated window
  DEBUG("SPECTRALREASSIGNMENT: Calculate time derivative window...");
  MatrixXR windowDeriv = _windowAlgo.window();
  for(int i = windowDeriv.cols() - 1; i > 0; i--){
    windowDeriv(0, i) = (windowDeriv(0, i) - windowDeriv(0, i - 1)) / timestep;
  }

  // TODO: Check what is the initial condition for the window
  // Should this be 0 or just the value it was originally * dt
  //windowDeriv(0, 0) = 0.0;
  _windowDerivAlgo.setWindow(windowDeriv);

  // Create the necessary buffers for the windowing
  _window.resize(1, _frameSize);
  _windowInteg.resize(1, _frameSize);
  _windowDeriv.resize(1, _frameSize);

  // Create the necessary buffers for the FFT
  _fft.resize(1, _fftSize);
  _fftAbs2.resize(1, _fftSize);
  _fftInteg.resize(1, _fftSize);
  _fftDeriv.resize(1, _fftSize);
  
  // Setup the algos
  _windowIntegAlgo.setup();
  _windowDerivAlgo.setup();
  _fftAlgo.setup();

  DEBUG("SPECTRALREASSIGNMENT: Finished set up...");
}


template<class F, class W, class ScalarW>
void SpectralReassignment::process(F frames, W* reassigned, W* fft){

  (*reassigned).resize(frames.rows(), _fftSize);
  (*fft).resize(frames.rows(), _fftSize);

  (*reassigned) = W::Zero((*reassigned).rows(), (*reassigned).cols());

  for (int i = 0; i < frames.rows(); i++){

    // Process the windowing
    _windowAlgo.process(frames, &_window);
    _windowIntegAlgo.process(frames, &_windowInteg);
    _windowDerivAlgo.process(frames, &_windowDeriv);

    // Process the FFT
    _fftAlgo.process(_window, &_fft);
    _fftAlgo.process(_windowInteg, &_fftInteg);
    _fftAlgo.process(_windowDeriv, &_fftDeriv);

    // Create the reassignment operations
    _fftAbs2 = (_fft.cwise().abs2()).template cast<Real>();
    
    DEBUG("SPECTRALREASSIGNMENT: Processing: creating the time reassignment operation...");    
    _reassignTime = -((_fftInteg.cwise() * _fft.conjugate()).cwise() / _fftAbs2.template cast<ScalarW>()).real().template cast<Real>();
    
    // TODO: Check the unity of the freq reassignment, it may need to be normalized by something
    DEBUG("SPECTRALREASSIGNMENT: Processing: creating the freq reassignment operation...");
    _reassignFreq = _freq + ((_fftDeriv.cwise() * _fft.conjugate()).cwise() / _fftAbs2.template cast<ScalarW>()).imag().template cast<Real>();
    
    // Reassign the spectrum values
    // TODO: put this into a function and do it right
    // will have to calculate and return all the reassigned values:
    // reassignedFrequency, reassignedTime:
    //      - are calculated using Flandrin's method using the 3 DFT
    // reassignedMagnitude, reassignedPhase: 
    //      - are calculated from reassigned freq and time and the original DFT
    //        (the magnitude and phase must then be put back 
    //         in the form of a complex in the reassigned frame)  
    DEBUG("SPECTRALREASSIGNMENT: Processing: reassigning...");
    DEBUG("SPECTRALREASSIGNMENT: Processing: reassigning _reassignFreq: " << _reassignFreq.rows() << ", " << _reassignFreq.cols());

    for(int j = 0; j < _reassignFreq.cols(); j++){

      if((int)round(_reassignFreq(i, j)) >= 0 && (int)round(_reassignFreq(i, j)) < (*reassigned).cols()) {

        (*reassigned)(i, (int)round(_reassignFreq(i, j))) += ((1.0 - (abs(_reassignFreq(i, j) - (int)round(_reassignFreq(i,j))))) * abs(_fft(i, (int)round(_reassignFreq(i,j)))));

      }
    }
  }
}

void SpectralReassignment::process(MatrixXC frames, MatrixXC* reassigned, MatrixXC* fft){
  process<MatrixXC, MatrixXC, Complex>(frames, reassigned, fft);
}

void SpectralReassignment::process(MatrixXR frames, MatrixXC* reassigned, MatrixXC* fft){
  process<MatrixXR, MatrixXC, Complex>(frames, reassigned, fft);
}

void SpectralReassignment::process(MatrixXC frames, MatrixXC* reassigned){
  process<MatrixXC, MatrixXC, Complex>(frames, reassigned, &_fft);
}

void SpectralReassignment::process(MatrixXR frames, MatrixXC* reassigned){
  process<MatrixXR, MatrixXC, Complex>(frames, reassigned, &_fft);
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
