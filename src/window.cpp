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

#include "debug.h"
#include "typedefs.h"

using namespace std;

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

Window::Window(int frameSize, Window::WindowType windowType) {
  DEBUG("WINDOW: Constructor frameSize: " << frameSize << ", windowType: " << windowType);

  _frameSize = frameSize;
  _windowType = windowType;

  DEBUG("WINDOW: Constructed");
}

Window::~Window(){
}

void Window::setup(){
  DEBUG("WINDOW: Setting up...");
  
  switch(_windowType){

  case RECTANGULAR:
    _window.set(MatrixXR::Ones(1, _frameSize));
    break;

  case HANN:
  case HANNING:
    _window.set(hann(_frameSize));
    break;

  case HAMMING:
    _window.set(hamming(_frameSize));
    break;

  case COSINE:
    _window.set(hamming(_frameSize));
    break;
    
  case BLACKMAN:
    _window.set(blackman(_frameSize));
    break;

  case NUTTALL:
    _window.set(nuttall(_frameSize));
    break;

  case BLACKMANHARRIS:
    _window.set(blackmanHarris(_frameSize));
    break;

  case BLACKMANNUTTALL:
    _window.set(blackmanNuttall(_frameSize));
    break;

  default:
    // Throw assertion unknown window type
    break;
  }
  
  DEBUG("WINDOW: Finished set up...");
}

MatrixXR Window::hann(int length){
  MatrixXR result(1, length);

  for(int i = 0; i < length; i++ ){
    result(0, i) = 0.5 * (1 - cos((2.0 * M_PI * (Real)i) / ((Real)length - 1.0)));
  }

  return result;
}

MatrixXR Window::hamming(int length){
  MatrixXR result(1, length);

  for(int i = 0; i < length; i++ ){
    result(0, i) = 0.53836 - 0.46164 * cos((2.0 * M_PI * (Real)i) / ((Real)length - 1.0));
  }

  return result;
}

MatrixXR Window::cosine(int length){
  MatrixXR result(1, length);

  for(int i = 0; i < length; i++ ){
    result(0, i) = sin((M_PI * (Real)i) / ((Real)length - 1.0));
  }

  return result;
}

MatrixXR Window::blackmanType(int length, Real a0, Real a1, Real a2, Real a3){
  MatrixXR result(1, length);

  Real pi_length_1 = M_PI / ((Real)length - 1.0);

  for(int i = 0; i < length; i++ ){
    result(0, i) = a0 \
                   - a1 * cos(2.0 * (Real)i * pi_length_1) \
                   + a2 * cos(4.0 * (Real)i * pi_length_1) \
                   + a3 * cos(6.0 * (Real)i * pi_length_1);
  }

  return result;
}

MatrixXR Window::blackman(int length){
  Real a0 = (1 - 0.16) / 2.0;
  Real a1 = 0.5;
  Real a2 = 0.16 / 2.0;
  Real a3 = 0.0;
  
  return blackmanType(length, a0, a1, a2, a3);
}


MatrixXR Window::nuttall(int length){
  Real a0 = 0.355768;
  Real a1 = 0.487396;
  Real a2 = 0.144232;
  Real a3 = 0.012604;
  
  return blackmanType(length, a0, a1, a2, a3);
}


MatrixXR Window::blackmanHarris(int length){
  Real a0 = 0.35875;
  Real a1 = 0.48829;
  Real a2 = 0.14128;
  Real a3 = 0.01168;
  
  return blackmanType(length, a0, a1, a2, a3);
}


MatrixXR Window::blackmanNuttall(int length){
  Real a0 = 0.3635819;
  Real a1 = 0.4891775;
  Real a2 = 0.1365995;
  Real a3 = 0.0106411;
  
  return blackmanType(length, a0, a1, a2, a3);
}

template<class F, class W>
void Window::process(F frames, W* windowedFrames){
  for (int i = 0; i < frames.rows(); i++){
    // Process and set
    (*windowedFrames).row(i) = frames.cwise() * _window;
  }
}

void Window::process(MatrixXC frames, MatrixXC* windowedFrames){
  process<MatrixXC, MatrixXC>(frames, windowedFrames);
}

void Window::process(MatrixXR frames, MatrixXC* windowedFrames){
  process<MatrixXR, MatrixXC>(frames, windowedFrames);
}

void Window::process(MatrixXR frames, MatrixXR* windowedFrames){
  process<MatrixXR, MatrixXR>(frames, windowedFrames);
}


void Window::reset(){
}

int Window::frameSize() const{
  return _frameSize;
}

Window::WindowType Window::windowType() const{
  return _windowType;
}

MatrixXR Window::window() const{
  return _window;
}
